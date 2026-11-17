"""APE-only check of the unchanged C-LPSD D (PSD) response to known tones.

Does not execute C, change its estimator, or invent an analytical Kaiser kernel.
Only consumes independent tone-only C outputs at two recorded peak amplitudes.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import h5py
import numpy as np


class KernelValidationError(ValueError):
    """Missing or inconsistent native C calibration evidence."""


def _sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as file:
        for block in iter(lambda: file.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def make_tone_fixture(path, *, frequency_hz, peak_amplitude_rad, sample_rate_hz,
                      duration_s, channel, phase_rad=0.0, chunk=262144):
    """Produce a deterministic single-channel radian cosine, not an axion signal."""
    path = Path(path).expanduser().resolve()
    values = (frequency_hz, peak_amplitude_rad, sample_rate_hz, duration_s, phase_rad)
    if channel not in ('X', 'Y') or not all(np.isfinite(values)):
        raise KernelValidationError('Finite tone parameters and channel X/Y required')
    if not (0 < frequency_hz < sample_rate_hz / 2 and peak_amplitude_rad > 0
            and sample_rate_hz > 0 and duration_s > 0 and chunk > 0):
        raise KernelValidationError('Invalid tone frequency, peak amplitude, duration or sample rate')
    n = int(round(sample_rate_hz * duration_s))
    if n < 4 or not np.isclose(n / sample_rate_hz, duration_s, rtol=0, atol=1e-10):
        raise KernelValidationError('Duration must contain an integer number of samples')
    if path.exists():
        raise KernelValidationError('Refusing to overwrite tone fixture')
    path.parent.mkdir(parents=True, exist_ok=True)
    try:
        with h5py.File(path, 'x') as h5:
            dataset = h5.create_dataset('rho_' + channel.lower(), shape=(n,), dtype='f8',
                                        chunks=(min(chunk, n),))
            for start in range(0, n, chunk):
                idx = np.arange(start, min(start + chunk, n), dtype=np.float64)
                dataset[start:start+len(idx)] = peak_amplitude_rad * np.cos(
                    2*np.pi*frequency_hz*idx/sample_rate_hz + phase_rad)
            for k, v in dict(frequency_hz=frequency_hz, peak_amplitude_rad=peak_amplitude_rad,
                             sample_rate_hz=sample_rate_hz, duration_s=duration_s,
                             phase_rad=phase_rad, channel=channel,
                             kind='SYNTHETIC_TONE_ONLY_NOT_AXION').items():
                h5.attrs[k] = v
    except BaseException:
        path.unlink(missing_ok=True)
        raise
    return {'file': str(path), 'sha256': _sha256(path), 'n_samples': n,
            'input_variance_rad2': peak_amplitude_rad**2/2}


def _selected_psd(path, *, native_frequency_hz, indices, columns):
    """Stream method-2 fDNbri (6 columns) or its lossless four-column combination."""
    target = set(indices.tolist())
    found, previous, count = {}, None, 0
    with Path(path).open('r', encoding='utf-8') as stream:
        for line_no, line in enumerate(stream, 1):
            if not line.strip() or line.lstrip().startswith('#'):
                continue
            fields = line.split()
            if len(fields) != columns:
                raise KernelValidationError(f'{path}:{line_no}: expected {columns} C columns')
            try:
                freq, psd, navg, nkernel = map(float, fields[:4])
            except ValueError as exc:
                raise KernelValidationError(f'{path}:{line_no}: nonnumeric C result') from exc
            if not all(np.isfinite((freq, psd, navg, nkernel))):
                raise KernelValidationError(f'{path}:{line_no}: nonfinite C result')
            if (freq <= 0 or psd < 0 or navg <= 0 or nkernel <= 0
                    or (previous is not None and freq <= previous)):
                raise KernelValidationError(f'{path}:{line_no}: invalid C frequency, PSD or metadata')
            if count >= len(native_frequency_hz) or freq != native_frequency_hz[count]:
                raise KernelValidationError(f'{path}:{line_no}: native C frequency grid mismatch')
            if count in target:
                found[count] = (psd, navg, nkernel)
            previous = freq
            count += 1
    if count != len(native_frequency_hz) or len(found) != len(indices):
        raise KernelValidationError(f'{path}: missing or extra native C bins')
    return np.array([found[i] for i in indices], dtype=np.float64)


def validate_c_kernel(*, native_frequency_hz, selected_indices, mode_frequencies_hz,
                      c_outputs, peak_amplitudes_rad, columns=6,
                      amplitude_rtol=0.02, active_floor_fraction=1e-6):
    """Measure K_D = D/(A_peak**2/2) and test quadratic scaling independently.

    c_outputs: {'X': [[low_file, high_file] per mode], 'Y': same}.
    Two C runs must use identical settings and phase except tone amplitude;
    this is an EXTERNAL provenance requirement, not established by C data alone.
    """
    native = np.asarray(native_frequency_hz, dtype=np.float64)
    idx = np.asarray(selected_indices)
    modes = np.asarray(mode_frequencies_hz, dtype=np.float64)
    amps = np.asarray(peak_amplitudes_rad, dtype=np.float64)
    if (native.ndim != 1 or native.size < 2 or not np.isfinite(native).all()
            or np.any(native <= 0) or np.any(np.diff(native) <= 0)):
        raise KernelValidationError('Invalid native C frequency grid')
    if (idx.ndim != 1 or idx.size == 0 or idx.dtype.kind not in 'iu'
            or np.any(idx < 0) or np.any(idx >= native.size) or np.any(np.diff(idx) <= 0)):
        raise KernelValidationError('Selected indices must be increasing native-grid integers')
    if (modes.ndim != 1 or not modes.size or not np.isfinite(modes).all()
            or np.any(modes <= 0) or np.any(np.diff(modes) <= 0)):
        raise KernelValidationError('Invalid sorted tone-mode frequencies')
    if (amps.shape != (2,) or not np.isfinite(amps).all() or amps[0] <= 0
            or amps[1] <= amps[0] or columns not in (4, 6)):
        raise KernelValidationError('Require two increasing positive peak amplitudes and 4/6-column format')
    if not (0 < amplitude_rtol < 1 and 0 < active_floor_fraction < 1):
        raise KernelValidationError('Invalid numerical validation tolerances')
    if set(c_outputs) != {'X', 'Y'} or any(len(c_outputs[ch]) != len(modes) for ch in ('X', 'Y')):
        raise KernelValidationError('Need independent X/Y paired C runs for every mode')
    kernel = np.empty((2, idx.size, modes.size), dtype=np.float64)
    diag = []
    for channel_idx, channel in enumerate(('X', 'Y')):
        for mode_idx, fmode in enumerate(modes):
            files = c_outputs[channel][mode_idx]
            if len(files) != 2 or files[0] == files[1]:
                raise KernelValidationError('Each tone mode needs distinct low/high C output paths')
            low = _selected_psd(files[0], native_frequency_hz=native,
                                indices=idx, columns=columns)
            high = _selected_psd(files[1], native_frequency_hz=native,
                                 indices=idx, columns=columns)
            if not np.array_equal(low[:, 1:], high[:, 1:]):
                raise KernelValidationError('C averaging/kernel metadata differ between tone amplitudes')
            k1 = low[:, 0] / (amps[0]**2 / 2)
            k2 = high[:, 0] / (amps[1]**2 / 2)
            maximum = max(float(k1.max()), float(k2.max()))
            if maximum <= 0:
                raise KernelValidationError('C response is zero for the selected tone bins')
            active = np.maximum(k1, k2) >= active_floor_fraction * maximum
            error = np.abs(k1[active] - k2[active]) / np.maximum(k1[active], k2[active])
            deviation = float(np.max(error))
            if deviation > amplitude_rtol:
                raise KernelValidationError(
                    f'{channel} tone {fmode:g} Hz fails C PSD amplitude-squared parity: {deviation:.3g}')
            # A second amplitude checks linearity ONLY; no invented correction or mean.
            kernel[channel_idx, :, mode_idx] = k1
            diag.append({'channel': channel, 'mode_hz': float(fmode),
                         'max_fractional_kernel_disagreement': deviation,
                         'active_bins_tested': int(np.count_nonzero(active)),
                         'input_files': [str(p) for p in files],
                         'input_sha256': [_sha256(p) for p in files]})
    return kernel, diag


def export_kernel(path, *, native_frequency_hz, selected_indices,
                  mode_frequencies_hz, c_outputs, peak_amplitudes_rad,
                  c_configuration_source, c_executable_sha256, **validation_options):
    """Save *provisional* native-C kernel in the exact Step-5 input array shape."""
    if not isinstance(c_configuration_source, str) or not c_configuration_source.strip():
        raise KernelValidationError('Native C command/configuration provenance required')
    if (not isinstance(c_executable_sha256, str) or len(c_executable_sha256) != 64
            or any(c not in '0123456789abcdefABCDEF' for c in c_executable_sha256)):
        raise KernelValidationError('Recorded C executable SHA-256 must be 64 hex digits')
    out = Path(path).expanduser().resolve()
    if out.suffix != '.npz' or out.exists():
        raise KernelValidationError('Refusing to overwrite existing/non-NPZ kernel')
    result, diagnostics = validate_c_kernel(
        native_frequency_hz=native_frequency_hz,
        selected_indices=selected_indices, mode_frequencies_hz=mode_frequencies_hz,
        c_outputs=c_outputs, peak_amplitudes_rad=peak_amplitudes_rad,
        **validation_options)
    metadata = {'status': 'PROVISIONAL_NATIVE_C_TONE_RESPONSE_NOT_PHYSICAL_APE_APPROVAL',
                'source_format': 'native C method-2 D = PSD rad^2/Hz',
                'kernel_units': 'Hz^-1 per unit time-domain rad^2 variance',
                'input_tone_peak_amplitudes_rad': list(map(float, peak_amplitudes_rad)),
                'c_configuration_source': c_configuration_source,
                'c_executable_sha256_declared_not_verified': c_executable_sha256,
                'c_provenance_check': 'C argv and input fixture evidence require separate verification',
                'diagnostics': diagnostics}
    out.parent.mkdir(parents=True, exist_ok=True)
    try:
        with out.open('xb') as file:
            np.savez(file, native_frequency_hz=np.asarray(native_frequency_hz),
                     output_bin_indices=np.asarray(selected_indices, dtype=np.int64),
                     halo_frequency_hz=np.asarray(mode_frequencies_hz),
                     native_c_psd_kernel_per_variance_hz=result,
                     metadata_json=json.dumps(metadata, sort_keys=True))
    except BaseException:
        out.unlink(missing_ok=True)
        raise
    return metadata


def load_campaign(path):
    """Load a bounded campaign; read native 5M-bin frequencies from Step-1 HDF5."""
    path = Path(path).expanduser().resolve()
    campaign = json.loads(path.read_text(encoding='utf-8'))
    if 'native_spectrum_hdf5' in campaign:
        if 'native_frequency_hz' in campaign:
            raise KernelValidationError('Choose native spectrum HDF5 or explicit frequencies, not both')
        source = Path(campaign.pop('native_spectrum_hdf5')).expanduser()
        if not source.is_absolute():
            source = (path.parent / source).resolve()
        with h5py.File(source, 'r') as h5:
            campaign['native_frequency_hz'] = h5['frequency'][:]
    return campaign


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--campaign-json', type=Path, required=True,
                        help='JSON with native_frequency_hz, selected_indices, mode_frequencies_hz, '
                             'c_outputs (X/Y x modes x [low,high]), peak_amplitudes_rad, '
                             'columns, c_configuration_source, c_executable_sha256')
    parser.add_argument('--output-npz', type=Path, required=True)
    options = parser.parse_args()
    campaign = load_campaign(options.campaign_json)
    print(json.dumps(export_kernel(options.output_npz, **campaign), indent=2))


if __name__ == '__main__':
    main()
