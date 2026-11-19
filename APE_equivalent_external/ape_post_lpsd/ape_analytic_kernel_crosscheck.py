#!/usr/bin/env python3
"""Compare an independently integrated Kaiser-window tone response to saved native C pilot outputs.

Diagnostic only. Does not edit native C or assert production physics/calibration.
C method 2's finite FFT interpolation and shared FFT approximation are NOT replicated.
"""
from pathlib import Path
import numpy as np
from numpy.polynomial.legendre import leggauss

BASE = Path.home() / 'Desktop/APE_Research'
FS, FMIN, FMAX, JDES, PSLL = 20000., 65., 10000., 5035956, 238.13
M = round(1/np.expm1(np.log(FMAX/FMIN)/(JDES-1)))


def response(f, f_tone):
    # Coefficients copied from original vendor LPSD/genwin.c:kaiser_alpha.
    x = PSLL/100
    alpha = ((0.0889732*x - 0.493285)*x + 4.71469)*x - 0.0821377
    L = np.floor(FS*M/f + 0.5).astype(np.int64)
    L -= L % 2
    f_eff = FS*M/L
    # C uses one Kaiser window per native bin, with length L.
    # Independent time-domain Fourier transform of this window, evaluated analytically.
    nodes, weights = leggauss(300)
    z = np.pi*alpha*np.sqrt(np.maximum(0., 1. - nodes**2))
    # Scale all Bessel values by I0(pi*alpha) to avoid large sums.
    # np.i0 arguments ~30 and are numerically safe.
    scaled = np.i0(z)/np.i0(np.pi*alpha)
    mean_w2 = 0.5*np.dot(weights, scaled**2)
    u = L*(f_tone-f_eff)/FS
    root2 = alpha**2-u**2
    kernel_ft = np.empty_like(u, dtype=np.float64)
    positive = root2 > 1e-12
    negative = root2 < -1e-12
    middle = ~(positive | negative)
    a = np.pi*np.sqrt(root2[positive])
    kernel_ft[positive] = np.sinh(a)/a/np.i0(np.pi*alpha)
    a = np.pi*np.sqrt(-root2[negative])
    kernel_ft[negative] = np.sin(a)/a/np.i0(np.pi*alpha)
    kernel_ft[middle] = 1./np.i0(np.pi*alpha)
    # For a cosine of peak amplitude A, variance A^2/2: PSD / variance = L/fs * |mean(W)|^2 / mean(W^2).
    return (L/FS)*kernel_ft**2/mean_w2


def main():
    all_good = True
    for label, directory, name in [
        ('5 kHz', 'APE_C_PSD_5kHz_Mac_pilot', 'mac_kernel_5kHz_DIAGNOSTIC_ONLY.npz'),
        ('8 kHz', 'APE_C_PSD_8kHz_Mac_pilot', 'mac_kernel_8kHz_DIAGNOSTIC_ONLY.npz')]:
        path = BASE/directory/name
        if not path.exists():
            print(f'{label}: MISSING {path}')
            all_good = False
            continue
        with np.load(path, allow_pickle=False) as z:
            freq = z['frequencies_hz'].astype(float)
            measured = z['K_D_per_variance'].astype(float)
        tone = float(freq[np.argmax(measured)])
        predicted = response(freq, tone)
        active = measured > measured.max()*1e-4
        relative = np.abs(predicted[active]-measured[active])/measured[active]
        peak = int(np.argmax(measured))
        print(f'\n{label}: synthetic tone at {tone:.9f} Hz')
        print(f'  C peak K_D:          {measured[peak]:.12g} Hz^-1')
        print(f'  analytic peak K_D:   {predicted[peak]:.12g} Hz^-1')
        print(f'  peak rel difference: {(predicted[peak]/measured[peak]-1)*100:+.5f}%')
        print(f'  active bins:         {active.sum()} of {len(freq)}')
        print(f'  active max rel diff: {relative.max()*100:.5f}%')
        print(f'  active median diff:  {np.median(relative)*100:.5f}%')
    print('\nSTATUS: INDEPENDENT ANALYTIC KAISER COMPARISON ONLY; not a full C kernel or APE physics validation')
    if not all_good: raise SystemExit(2)

if __name__=='__main__': main()