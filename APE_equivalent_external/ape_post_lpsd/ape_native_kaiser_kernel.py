#!/usr/bin/env python3
"""Native-grid analytic Kaiser PSD tone response for APE C-LPSD.

Source model: original vendor genwin.c Kaiser and lpsd.c method-2 window
normalization. This computes a deterministic finite-window expectation, NOT
an exact reproduction of the method-2 truncated FFT interpolation, a physical
halo, or a validated coupling sensitivity. NO vendor code is changed.

Inputs: actual native C-LPSD frequency dataset; one or several explicitly
requested tone frequencies. Computes K_D [Hz^-1 / unit time-domain variance]
for a limited set of native output bins around each tone; never exports a full
5-million-bin physical response or assigns axion parameters.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
import h5py
from numpy.polynomial.legendre import leggauss


def kaiser_alpha(psll):
    x = float(psll) / 100.0
    return ((0.0889732 * x - 0.493285) * x + 4.71469) * x - 0.0821377


def native_lengths(frequencies, fs, fmin, fmax, jdes):
    m = round(1.0 / np.expm1(np.log(fmax / fmin) / (jdes - 1)))
    L = np.floor(fs * m / np.asarray(frequencies, dtype=np.float64) + 0.5).astype(np.int64)
    L -= L % 2
    if np.any(L <= 0):
        raise ValueError('Nonpositive native window length')
    return m, L


def model_kernel(output_frequencies, mode_frequencies, *, fs, fmin, fmax, jdes, psll=238.13):
    """Return K [native-output bins, supplied input modes], no bin renormalization."""
    f = np.asarray(output_frequencies, np.float64)
    modes = np.asarray(mode_frequencies, np.float64)
    if f.ndim != 1 or modes.ndim != 1 or not len(f) or not len(modes):
        raise ValueError('Frequencies and modes must be nonempty 1D arrays')
    if not (np.isfinite(f).all() and np.isfinite(modes).all() and
            (f > 0).all() and (modes > 0).all()):
        raise ValueError('Frequencies must be finite and positive')
    if not 25 <= psll <= 250:
        raise ValueError('Unsupported Kaiser PSLL')
    m, L = native_lengths(f, fs, fmin, fmax, jdes)
    eff = fs * m / L
    alpha = kaiser_alpha(psll)
    nodes, weights = leggauss(300)
    denom = np.i0(np.pi * alpha)
    w = np.i0(np.pi * alpha * np.sqrt(np.maximum(0, 1 - nodes**2))) / denom
    mean_w2 = 0.5 * np.dot(weights, w * w)
    u = (L[:, None] / fs) * (modes[None, :] - eff[:, None])
    d = alpha * alpha - u * u
    ft = np.empty_like(d)
    pos = d > 1e-12
    neg = d < -1e-12
    mid = ~(pos | neg)
    a = np.pi * np.sqrt(d[pos]); ft[pos] = np.sinh(a) / a / denom
    a = np.pi * np.sqrt(-d[neg]); ft[neg] = np.sin(a) / a / denom
    ft[mid] = 1 / denom
    return (L[:, None] / fs) * (ft * ft) / mean_w2


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--spectrum', type=Path, required=True, help='Actual native C-LPSD HDF5')
    p.add_argument('--tone-hz', type=float, action='append', required=True,
                   help='Tone frequency; repeat to provide multiple modes')
    p.add_argument('--half-width-bins', type=int, default=40)
    p.add_argument('--psll', type=float, default=238.13, help='PROVISIONAL until production PSLL verified')
    p.add_argument('--output', type=Path, required=True, help='New diagnostic .npz (no overwrite)')
    a = p.parse_args()
    if a.half_width_bins < 1 or a.half_width_bins > 5000:
        p.error('Choose half-width-bins from 1 to 5000')
    if a.output.exists() or a.output.suffix != '.npz':
        p.error('Output must be a new .npz file')
    modes = np.array(a.tone_hz, dtype=np.float64)
    with h5py.File(a.spectrum, 'r') as h:
        ds = h['frequency']
        n = len(ds)
        attrs = h['logPSD'].attrs
        fs, fmin, fmax, jdes = (float(attrs['fs']), float(attrs['fmin']),
                               float(attrs['fmax']), int(attrs['Jdes']))
        if n != jdes:
            raise ValueError(f'Frequency grid length {n} != Jdes {jdes}')
        # The C native grid f_j = fmin exp(j ln(fmax/fmin)/(Jdes-1)).
        # Find native index algebraically, then verify using real HDF5 data.
        indices = []
        for tone in modes:
            j = int(round((jdes-1)*np.log(tone/fmin)/np.log(fmax/fmin)))
            if not a.half_width_bins <= j < n-a.half_width_bins:
                raise ValueError(f'Tone {tone:g} too close to grid boundary')
            indices.extend(range(j-a.half_width_bins, j+a.half_width_bins+1))
        indices = np.unique(np.asarray(indices, dtype=np.int64))
        # Read only the desired native windows, via contiguous slices.
        frequencies = np.empty(len(indices), dtype=np.float64)
        breaks = np.r_[0, 1+np.flatnonzero(np.diff(indices) != 1), len(indices)]
        for left, right in zip(breaks[:-1], breaks[1:]):
            frequencies[left:right] = ds[int(indices[left]):int(indices[right-1])+1]
        if not np.allclose(frequencies, fmin*np.exp(indices*np.log(fmax/fmin)/(jdes-1)),
                           rtol=1e-8, atol=1e-7):
            raise ValueError('Actual HDF5 frequency grid differs from original C native grid')
    kernel = model_kernel(frequencies, modes, fs=fs, fmin=fmin, fmax=fmax,
                          jdes=jdes, psll=a.psll)
    metadata = dict(status='DIAGNOSTIC_ANALYTIC_KAISER_NOT_PHYSICAL_AXION_RESPONSE',
                    source_spectrum=str(a.spectrum.resolve()),
                    psll_db=a.psll, psll_production_verified=False,
                    estimator='window-integral approximation to original C-LPSD method 2',
                    output_units='Hz^-1 per unit input time-domain variance',
                    shape='K_D[output_bin, input_tone]',
                    no_rebinning=True, no_normalization_to_unit_bin_integral=True)
    a.output.parent.mkdir(parents=True, exist_ok=True)
    np.savez(a.output, native_frequency_hz=frequencies,
             native_bin_indices=indices, input_tone_hz=modes,
             K_D_per_variance_hz=kernel,
             metadata_json=np.array(json.dumps(metadata, sort_keys=True)))
    print(json.dumps({'status':metadata['status'],'native_bins':int(len(indices)),
                      'input_modes':modes.tolist(), 'peak_k_hz_inv':kernel.max(axis=0).tolist(),
                      'output':str(a.output.resolve())}, indent=2))

if __name__ == '__main__':
    main()