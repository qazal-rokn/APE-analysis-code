# Paper 2026 to APE signal-model audit (28 October 2026)

**Scope:** source-based, read-only scientific audit. No C, vendor, calibration, or APE scientific code was changed. No new synthetic injection was generated.

**Source:** Göttel and Raymond, *Analytical kernels for efficient constant Q transforms in dark matter searches with LIGO*, Scientific Reports 16, 15364 (2026), https://doi.org/10.1038/s41598-025-33428-2. Compared with unchanged original `stats.py`, `max_lkl_executable.py`, `fit_background.py`, `calc_upper_limit.py`, and the APE development adapters.

## Findings supported by the paper

1. The **scalar** field and LIGO strain/coupling transfer (paper Eqs. 1–4) are not physically interchangeable with APE's pseudoscalar axion-photon QWP polarimetry response.
2. The paper assumes a virialized, isotropic/homogeneous Maxwell–Boltzmann velocity distribution, escape-speed cutoff, and solar-system and Earth orbital effects on broadening (Eq. 5 and text). It does not provide a validated time-domain APE axion X/Y injection generator.
3. Paper Eqs. 7–10 specify `Q=f/delta_f`, logarithmic bin spacing, and `N_j = fs/(f_{j+1}-f_j)`. For f=5000 Hz, Q=1e6 implies delta_f=0.005 Hz, and approximately 200 seconds per coherent C window.
4. Paper Eq. 11 computes window-normalized averaged spectral density with 20% overlap. The APE PS-enabled C output is **direct PS in rad^2**, not PSD in rad^2/Hz; units and native kernel need explicit handling.
5. Paper background fitting uses a local log-log spline and skew-normal *log-power* residuals, rejecting poor fits outside 0.5–10. The unchanged vendor background implementation passed numerical parity and selected native-C tests.
6. Paper Eqs. 16–17 and original vendor `stats.log_likelihood` use `log(observed PS)-log(background + mu * peak_norm * peak_shape)` with skew-normal residual likelihood. The paper explicitly **neglects dark-matter signal stochasticity in this likelihood because the chosen integration window matches coherence time**. It does **not** authorize invented independent random signal-power draws for every bin. The paper also has a LIGO-specific calibration nuisance term, not automatically the correct APE nuisance term.
7. Paper Eqs. 18–19 define one-sided `q0` and upper-limit `q_mu`; the original vendor expected limit is an asymptotic method. The tested APE adapters provide software parity, not APE global significance, observed limits or validated 95% coverage.

## Key mismatch: previously successful 5 kHz reduced-grid pilot

The software-integration pilot used `fs=20000 Hz`, `duration=200 s`, `fmin=4995 Hz`, `fmax=5008 Hz`, and `J=1301`. On that **actual grid**:

```text
Q_effective = 1/[exp(ln(5008/4995)/1300)-1] = 500149.2184
delta_f(5000 Hz) ≈ 0.00999701652 Hz
native C integration window at 5000 Hz ≈ 100.02984 s
```

**Verdict:** authentic compiled-C end-to-end execution, background, vendor MLE, local candidates and UL software integration **PASS** for the pilot; exact paper-coherence matching (`Q=1e6`) **NOT PASSED by this pilot**. Do not convert previous software acceptance into physical halo validation.

## Source-derived next Q=10^6 grid (not yet executed)

Use 1301 bins, center j=650 at 5000 Hz, `s=ln(1+1e-6)`:

```text
f_j = 5000 * exp((j-650)*s) Hz
fmin = 4996.751057645126 Hz
fmax = 5003.251054852836 Hz
Q_effective = 1,000,000
center delta_f = 0.005 Hz
center N_j at fs=20000 Hz = 4,000,000 samples = 200 s
lowest-frequency largest N_j ≈ 4,002,601 samples = 200.13005 s
```

A **201 s** recording supplies enough samples to contain at least one complete C window at every bin, subject to confirming actual native averaging/overlap settings. This is a proposed, resource-bounded C **coherent-window sanity check**, not a claim of Maxwellian-halo or physical APE validation; performance has not yet been measured.

## Remaining scientific contracts

- Derive/validate an APE-specific pseudoscalar `g_aγγ -> rho_x/rho_y` amplitude and phase response, including peak/RMS convention, joint X/Y dependence, QWP/cavity transfer, and measured native-C window/PS response.
- Source the exact intended halo model and field realization prescription; do not invent per-bin exponential signal-power fluctuations. Distinguish scalar-LIGO from pseudoscalar-APE physics.
- Verify single-bin versus multi-bin template validity, covariance between overlapping C bins, local/global candidate statistics, injection-strength ensemble performance and APE confidence coverage.
- Preserve `main`, raw/frozen data, compiled original C, original vendor and prior results read-only.

## Audit status

**Paper-to-code assumption inspection: complete.** **Physical APE signal-injection contract: not yet approved or validated.** No code changed by this audit.
