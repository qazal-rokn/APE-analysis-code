"""Lossless APE C-LPSD-to-vendor IO handoff; no statistical model changes.

Input: independent four-column combined spectra: f, PSD, N_avg, N_kernel.
Output: HDF5 accepted by the unchanged vendor LPSDOutput class.

The frequency grid and PSD are preserved as float64 arrays; only logPSD,
which the vendor reader also computes, is additionally stored.
"""
from __future__ import annotations

import argparse
from itertools import zip_longest
import hashlib
import json
import math
from pathlib import Path

import h5py
import numpy as np


class HandoffError(ValueError):
    """The requested handoff violates the recorded APE input contract."""


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open('rb') as source:
        for block in iter(lambda: source.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def _rows(path: Path):
    with path.open('r', encoding='utf-8') as stream:
        for line_number, line in enumerate(stream, 1):
            if not line.strip() or line.lstrip().startswith('#'):
                continue
            columns = line.split()
            if len(columns) != 4:
                raise HandoffError(f'{path}:{line_number}: expected exactly four columns (f, PSD, N_avg, N_kernel)')
            try:
                f, psd, n_avg, n_kernel = (float(item) for item in columns)
            except ValueError as exc:
                raise HandoffError(f'{path}:{line_number}: nonnumeric C-LPSD row') from exc
            if not all(math.isfinite(v) for v in (f, psd, n_avg, n_kernel)):
                raise HandoffError(f'{path}:{line_number}: nonfinite C-LPSD row')
            if f <= 0 or psd <= 0 or n_avg <= 0 or n_kernel < 0:
                raise HandoffError(f'{path}:{line_number}: invalid frequency, PSD or averaging metadata')
            yield f, psd, n_avg, n_kernel


def _new_hdf5(path: Path, *, chunk_size: int, channel: str, fs: float, t0: int, t1: int):
    h5 = h5py.File(path, 'x')
    try:
        h5.attrs['channel'] = channel
        h5.attrs['source_format'] = 'combined_ape_c_method2_f_PSD_N_avg_N_kernel'
        h5.attrs['spectrum_unit'] = 'rad^2/Hz'
        h5.attrs['sample_rate_hz'] = fs
        h5.attrs['segment_t0'] = t0
        h5.attrs['segment_t1'] = t1
        for key in ('frequency', 'PSD', 'logPSD', 'N_avg', 'N_kernel'):
            h5.create_dataset(key, shape=(0,), maxshape=(None,), chunks=(chunk_size,), dtype='f8')
    except BaseException:
        h5.close()
        raise
    return h5


def _append(h5: h5py.File, rows: list[tuple[float, float, float, float]]) -> None:
    if not rows:
        return
    values = np.asarray(rows, dtype=np.float64)
    n = len(values)
    arrays = {
        'frequency': values[:, 0],
        'PSD': values[:, 1],
        'logPSD': np.log(values[:, 1]),
        'N_avg': values[:, 2],
        'N_kernel': values[:, 3],
    }
    for name, array in arrays.items():
        ds = h5[name]
        start = ds.shape[0]
        ds.resize((start + n,))
        ds[start:start + n] = array


def prepare_pair(x_path: str | Path, y_path: str | Path, output_dir: str | Path, *,
                 fs: float, t0: int, t1: int, expected_bins: int | None = None,
                 chunk_size: int = 16384) -> dict:
    """Create two original-reader-compatible HDF5 products, without modifying C data.

    Never combine X/Y, never rebin or rescale PSD, and never change the vendor
    source. Both files must have the exact same C frequency grid.
    """
    x_path, y_path, output_dir = map(lambda p: Path(p).expanduser().resolve(), (x_path, y_path, output_dir))
    if x_path == y_path:
        raise HandoffError('X/Y must be different source files')
    if not all(p.is_file() for p in (x_path, y_path)):
        raise HandoffError('Both C-LPSD input files must exist')
    if not math.isfinite(fs) or fs <= 0:
        raise HandoffError('Sample rate must be positive and finite')
    if not isinstance(t0, int) or not isinstance(t1, int) or t0 < 0 or t1 <= t0:
        raise HandoffError('Provide the recorded positive integer segment start/end identifiers (t0 < t1)')
    if expected_bins is not None and expected_bins <= 0:
        raise HandoffError('Expected bin count must be positive')
    if chunk_size <= 0:
        raise HandoffError('Chunk size must be positive')
    if output_dir == x_path or output_dir == y_path or output_dir in x_path.parents or output_dir in y_path.parents:
        raise HandoffError('Output directory must not contain the source files')

    names = {channel: f'APE_{t0}_{t1}_{channel}.h5' for channel in ('X', 'Y')}
    output = {channel: output_dir / name for channel, name in names.items()}
    report_file = output_dir / f'APE_{t0}_{t1}_handoff.json'
    staging = {channel: path.with_suffix(path.suffix + '.partial') for channel, path in output.items()}
    for path in (*output.values(), *staging.values(), report_file):
        if path.exists():
            raise HandoffError(f'Refusing to overwrite existing destination: {path}')

    output_dir.mkdir(parents=True, exist_ok=True)
    h5s: dict[str, h5py.File] = {}
    completed = []
    count = 0
    fmin = fmax = previous = None
    xbuf: list[tuple[float, float, float, float]] = []
    ybuf: list[tuple[float, float, float, float]] = []
    try:
        for channel in ('X', 'Y'):
            h5s[channel] = _new_hdf5(staging[channel], chunk_size=chunk_size,
                                     channel=channel, fs=fs, t0=t0, t1=t1)
        for x_row, y_row in zip_longest(_rows(x_path), _rows(y_path)):
            if x_row is None or y_row is None:
                raise HandoffError(f'X/Y frequency bin count differs after {count} rows')
            fx, fy = x_row[0], y_row[0]
            if fx != fy:
                raise HandoffError(f'X/Y native frequency mismatch at row {count + 1}: {fx} != {fy}')
            if previous is not None and fx <= previous:
                raise HandoffError(f'Non-increasing/duplicate native frequency at row {count + 1}')
            if fx > fs / 2.0:
                raise HandoffError(f'Frequency {fx} exceeds the declared Nyquist limit {fs/2.0}')
            if fmin is None:
                fmin = fx
            fmax = previous = fx
            xbuf.append(x_row)
            ybuf.append(y_row)
            count += 1
            if len(xbuf) == chunk_size:
                _append(h5s['X'], xbuf)
                _append(h5s['Y'], ybuf)
                xbuf.clear()
                ybuf.clear()
        if count < 2:
            raise HandoffError('Need at least two native frequency bins')
        if expected_bins is not None and count != expected_bins:
            raise HandoffError(f'Expected {expected_bins} frequency bins; found {count}')
        _append(h5s['X'], xbuf)
        _append(h5s['Y'], ybuf)
        for channel in ('X', 'Y'):
            group = h5s[channel]['logPSD'].attrs
            # These are exactly the names read by the unchanged vendor LPSDOutput.
            group['fmin'] = fmin
            group['fmax'] = fmax
            group['fs'] = fs
            group['Jdes'] = count
            h5s[channel].attrs['source_sha256'] = _sha256(x_path if channel == 'X' else y_path)
            h5s[channel].attrs['native_frequency_grid'] = True
            h5s[channel].flush()
            h5s[channel].close()
        h5s.clear()

        summary = {
            'status': 'source-preserving APE-to-vendor LPSD IO handoff; no fitter executed',
            'input_columns': ['frequency', 'PSD', 'N_avg', 'N_kernel'],
            'vendor_datasets': ['frequency', 'logPSD'],
            'extra_preserved_datasets': ['PSD', 'N_avg', 'N_kernel'],
            'native_frequency_grid_exactly_shared': True,
            'independent_X_Y': True,
            'bin_count_per_quadrature': count,
            'fmin_hz': fmin,
            'fmax_hz': fmax,
            'sample_rate_hz': fs,
            't0': t0,
            't1': t1,
            'source': {'X': str(x_path), 'Y': str(y_path)},
            'source_sha256': {'X': _sha256(x_path), 'Y': _sha256(y_path)},
            'output': {channel: str(path) for channel, path in output.items()},
        }
        for channel in ('X', 'Y'):
            staging[channel].rename(output[channel])
            completed.append(output[channel])
        with report_file.open('x', encoding='utf-8') as stream:
            json.dump(summary, stream, indent=2)
            stream.write('\n')
        return summary
    except BaseException:
        for h5 in h5s.values():
            h5.close()
        for path in (*staging.values(), *completed):
            path.unlink(missing_ok=True)
        raise


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--x-psd', required=True, help='Four-column rho_x combined C-LPSD spectrum')
    p.add_argument('--y-psd', required=True, help='Four-column rho_y combined C-LPSD spectrum')
    p.add_argument('--output-dir', required=True, help='Non-overwriting output directory')
    p.add_argument('--fs', required=True, type=float, help='Original sampling frequency in Hz (APE: usually 20000)')
    p.add_argument('--t0', required=True, type=int, help='Actual source-segment start identifier (integer time)')
    p.add_argument('--t1', required=True, type=int, help='Actual source-segment stop identifier (integer time)')
    p.add_argument('--expected-bins', type=int, default=None, help='Historical APE Constant-Q: 5035956')
    p.add_argument('--chunk-size', type=int, default=16384)
    a = p.parse_args()
    print(json.dumps(prepare_pair(a.x_psd, a.y_psd, a.output_dir, fs=a.fs, t0=a.t0, t1=a.t1,
                                  expected_bins=a.expected_bins, chunk_size=a.chunk_size), indent=2))


if __name__ == '__main__':
    main()
