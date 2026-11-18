#!/usr/bin/env python3
"""APE four-column C-LPSD -> HDF5 using the UNMODIFIED external LPSDOutput.

The adapter ONLY adds the vendor text header and invokes the original reader
and original LPSDData.save_to_hdf5. Thus the original log-grid reconstruction
is retained, deliberately. It does not modify likelihood or background code.

Requires the existing APE_equivalent_external/scalar_dark_matter_lpsd package.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
from itertools import zip_longest

import h5py
import numpy as np

# When installed alongside ape_post_lpsd scripts, locate the original package.
_VENDOR_ROOT = Path(__file__).resolve().parents[1]
if (_VENDOR_ROOT / "scalar_dark_matter_lpsd").is_dir():
    sys.path.insert(0, str(_VENDOR_ROOT))

from scalar_dark_matter_lpsd.LPSDIO import LPSDOutput  # ORIGINAL, unmodified


class VendorConversionError(ValueError):
    pass


def _scan_ape_file(path: Path) -> tuple[int, float, float, str]:
    """Inspect APE text input without replacing any scientific values."""
    h = hashlib.sha256()
    count, f_first, f_last = 0, None, None
    with path.open('rb') as inp:
        for n, raw in enumerate(inp, 1):
            h.update(raw)
            text = raw.decode('utf-8').strip()
            if not text or text.startswith('#'):
                continue
            cols = text.split()
            if len(cols) != 4:
                raise VendorConversionError(f'{path}:{n}: expected f PSD N_avg N_kernel')
            try:
                f, psd, n_avg, n_kernel = map(float, cols)
            except ValueError as exc:
                raise VendorConversionError(f'{path}:{n}: nonnumeric row') from exc
            if not all(map(math.isfinite, (f, psd, n_avg, n_kernel))):
                raise VendorConversionError(f'{path}:{n}: nonfinite row')
            if f <= 0 or psd <= 0 or n_avg <= 0 or n_kernel < 0:
                raise VendorConversionError(f'{path}:{n}: invalid spectral values')
            if f_last is not None and f <= f_last:
                raise VendorConversionError(f'{path}:{n}: nonincreasing frequency')
            if f_first is None:
                f_first = f
            f_last = f
            count += 1
    if count < 2:
        raise VendorConversionError(f'{path}: requires at least two spectral bins')
    return count, f_first, f_last, h.hexdigest()


def _stage_vendor_text(source: Path, destination: Path, fmin: float,
                       fmax: float, fs: float, bins: int) -> None:
    # Preserve all four numerical columns; add only the vendor-required header.
    # The original text parser also requires tab-separated numeric rows.
    header = f'# Command line: lpsd -s {fmin:.17g} -t {fmax:.17g} -f {fs:.17g} -J {bins}\n'
    with destination.open('w', encoding='utf-8') as target:
        target.write(header)
        with source.open('r', encoding='utf-8') as inp:
            for line in inp:
                if line.strip() and not line.lstrip().startswith('#'):
                    target.write('\t'.join(line.split()) + '\n')


def _compare_grid(source: Path, reconstructed: np.ndarray) -> dict:
    """Diagnostic only: record vendor-reconstructed grid vs C frequencies."""
    i = 0
    max_abs = max_rel = 0.0
    mismatches = 0
    with source.open('r', encoding='utf-8') as inp:
        for line in inp:
            if not line.strip() or line.lstrip().startswith('#'):
                continue
            f = float(line.split()[0])
            d = abs(float(reconstructed[i]) - f)
            max_abs = max(max_abs, d)
            max_rel = max(max_rel, d / f)
            if d > 1e-10 * f:
                mismatches += 1
            i += 1
    return {'max_absolute_difference_hz': max_abs,
            'max_relative_difference': max_rel,
            'bins_differing_above_1e-10_relative': mismatches}


def _check_matching_native_grids(x: Path, y: Path) -> None:
    """Check only; do not alter the original reader's frequency reconstruction."""
    def frequencies(path):
        with path.open('r', encoding='utf-8') as file:
            for line in file:
                if line.strip() and not line.lstrip().startswith('#'):
                    yield float(line.split()[0])
    for i, (fx, fy) in enumerate(zip_longest(frequencies(x), frequencies(y)), 1):
        if fx != fy:
            raise VendorConversionError(f'X/Y native frequency mismatch at row {i}')


def convert_pair_original(x_psd: str | Path, y_psd: str | Path,
                          output_dir: str | Path, *, fs: float, t0: int, t1: int,
                          expected_bins: int | None = None) -> dict:
    """Use original vendor LPSDOutput and save_to_hdf5 for each APE quadrature."""
    x, y, out = (Path(p).expanduser().resolve() for p in (x_psd, y_psd, output_dir))
    if x == y or not x.is_file() or not y.is_file():
        raise VendorConversionError('X and Y must be distinct existing files')
    if not math.isfinite(fs) or fs <= 0 or not 0 <= t0 < t1:
        raise VendorConversionError('Invalid sampling rate or segment identifiers')
    a, b = _scan_ape_file(x), _scan_ape_file(y)
    if a[:3] != b[:3]:
        raise VendorConversionError('X/Y bin counts or endpoint frequencies differ')
    _check_matching_native_grids(x, y)
    if expected_bins is not None and a[0] != expected_bins:
        raise VendorConversionError(f'Expected {expected_bins} bins, found {a[0]}')
    if a[2] > fs / 2:
        raise VendorConversionError('Frequency exceeds Nyquist limit')
    # Original conversion retains only frequency and logPSD: no extra
    # normalization, native-grid override, or auxiliary HDF5 datasets.
    paths = {q: out / f'APE_{t0}_{t1}_{q}.h5' for q in ('X', 'Y')}
    report_path = out / f'APE_{t0}_{t1}_original_converter_report.json'
    if any(p.exists() for p in (*paths.values(), report_path)):
        raise FileExistsError('Refusing to overwrite an existing conversion output')
    out.mkdir(parents=True, exist_ok=True)
    report = {
        'converter': 'original_scalar_dark_matter_lpsd.LPSDIO.LPSDOutput',
        'frequency_rule': 'original vendor logspace reconstruction (NOT native-grid override)',
        'n_bins': a[0], 'fs_hz': fs, 't0': t0, 't1': t1,
        'inputs': {'X': str(x), 'Y': str(y)},
        'source_sha256': {'X': a[3], 'Y': b[3]},
        'outputs': {q: str(p) for q, p in paths.items()},
        'grid_comparison': {},
    }
    created = []
    try:
        with TemporaryDirectory(prefix='ape_vendor_original_') as temp:
            for q, src in (('X', x), ('Y', y)):
                staged = Path(temp) / f'vendor_input_{q}.lpsd'
                _stage_vendor_text(src, staged, a[1], a[2], fs, a[0])
                original = LPSDOutput(str(staged))
                report['grid_comparison'][q] = _compare_grid(src, original.freq)
                # Original code includes optional absent flags as None. h5py
                # cannot serialize None; omit ONLY those absent attributes.
                original.kwargs = {key: value for key, value in original.kwargs.items()
                                   if value is not None}
                created.append(paths[q])
                original.save_to_hdf5(str(paths[q]))  # ORIGINAL method
                with h5py.File(paths[q], 'r') as check:
                    if not np.array_equal(check['frequency'][:], original.freq):
                        raise RuntimeError('Vendor frequency round-trip failed')
                    if not np.array_equal(check['logPSD'][:], original.logPSD):
                        raise RuntimeError('Vendor PSD round-trip failed')
        with report_path.open('x', encoding='utf-8') as handle:
            json.dump(report, handle, indent=2)
            handle.write('\n')
        return report
    except BaseException:
        for file in created:
            file.unlink(missing_ok=True)
        report_path.unlink(missing_ok=True)
        raise


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--x-psd', required=True)
    p.add_argument('--y-psd', required=True)
    p.add_argument('--output-dir', required=True)
    p.add_argument('--fs', required=True, type=float)
    p.add_argument('--t0', required=True, type=int)
    p.add_argument('--t1', required=True, type=int)
    p.add_argument('--expected-bins', type=int, default=None)
    args = p.parse_args()
    result = convert_pair_original(args.x_psd, args.y_psd, args.output_dir,
                                   fs=args.fs, t0=args.t0, t1=args.t1,
                                   expected_bins=args.expected_bins)
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
