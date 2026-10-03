# Low-Frequency Full-Record FFT Extension

## Scope

This independent package is a future low-frequency PSD method. It accepts one
already validated, uninterrupted, ordered voltage sequence for **one** channel
at a time. It does not perform data loading, cleaning, calibration, segment
selection, gap filling, interpolation, resampling, reordering, or X/Y
combination. Its method identity is always `fft_full_record`; the compiled
estimator remains `c_lpsd`.

```text
validated pem1_x voltage ──> FFT extension ──> FFT PSD/ASD product
                         └─> C-LPSD backend ─> C-LPSD PSD/ASD product
validated pem1_Y voltage ──> FFT extension ──> FFT PSD/ASD product
                         └─> C-LPSD backend ─> C-LPSD PSD/ASD product

method-specific numerical products ──> optional display-only combined plot
```

The optional combined display never overwrites, concatenates, interpolates, or
otherwise merges the underlying FFT and C-LPSD numerical products.

## Planned FFT convention

The future kernel performs exactly one full-record `rfft` on a single voltage
sequence. It returns native FFT bins and a one-sided density PSD:

```text
Sxx[k] = |RFFT{x[n] w[n]}|^2 / (fs sum_n w[n]^2)
```

Interior positive-frequency bins are doubled; DC and (for even record length)
the Nyquist bin are not. Thus the intended units are `V^2/Hz` and
`V/sqrt(Hz)`. It is not Welch averaging, does not use shorter subwindows, and
does not claim to be identical to C-LPSD.

## C-window inspection and status

The unchanged C source was inspected. `genwin.c` maps PSLL to a C Kaiser
`alpha` polynomial and samples
`I0(pi alpha sqrt(1-z^2))/I0(pi alpha)` with `z=2i/(N-1)-1`. Its LPSD DFT
normalization uses the actual `sum(w^2)` (see `lpsd.c`, `getDFT2`:
`2/(fs sum(w^2))`). This establishes the C convention but does **not** prove
that NumPy's Kaiser implementation or its numeric Bessel evaluation is
identical. Therefore `numpy_kaiser_to_be_validated` is explicitly labelled
*to be validated*, while the default is an explicit rectangular window.

## Ranges and preflight

- **Requested display range:** future figure limits only.
- **FFT-calculable range:** native full-record bins from `fs/N` through Nyquist;
  configured FFT bounds select/report bins but do not alter the transform.
- **C-LPSD-calculable range:** range in the immutable C-native manifest.
- **Handoff:** a configured comparison/display boundary that must lie in both
  declared ranges and must later be checked against the actual C manifest.

Before any allocation, `preflight_full_record_fft` calculates a conservative
working-set estimate (input, detrended and windowed work arrays, window,
complex rFFT, PSD, and frequency arrays) and refuses if it exceeds the
configured memory limit. The FFT extension has no authority to relax that
limit or alter the record.

Use `configs/fft_extension.example.yaml` as the starting configuration. No
execution command is enabled or run as part of this implementation stage.

## Future no-data preflight

`tools/preflight_fft_extension.py --help` documents a future metadata-only
preflight. It requires caller-provided validated `--sample-count` and
`--sample-rate-hz`, reads no sequence values, and reports native bin width,
FFT-calculable range, and conservative expected memory before any allocation.
It must pass before a future FFT execution is considered.
