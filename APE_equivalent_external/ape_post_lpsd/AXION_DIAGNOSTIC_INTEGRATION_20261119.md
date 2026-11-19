# APE axion kernel and likelihood: 19 November 2026 integration

This integrates the **actual previously executed 5/8 kHz engineering scripts**
into the existing APE adapter package, without modifying upstream
Scalar-Dark-Matter-LPSD C code or Python likelihood, background fitting,
optimization, candidate selection, or upper-limit functions.

**Status: DIAGNOSTIC ONLY. Not an approved full-band axion search or
physical coupling upper limit.** Source-specific PSLL, X/Y optical projections,
phase response, statistical coverage, and halo assumptions remain open.

## Files

- `ape_native_kaiser_kernel.py`: analytic Kaiser one-sided D PSD
  response per unit input time-domain variance, on exact native output bins.
- `ape_axion_halo_kernel_diagnostic.py`: reproducible fixed-velocity
  benchmark halo and convolution with native analytic C kernel.
- `ape_qwp_halo_psd_preview.py`: reuses unchanged existing QWP physics
  function and records measured X/Y PSD alongside *total* template.
- `ape_provisional_qwp_multibin_vendor.py`: unchanged original vendor
  optimizer on 81 overlapping bins, separate X/Y counterfactual projections.
- `ape_signal_weighted_covariance_diagnostic.py`: empirical lag correlations,
  illustrative signal-weighted standard deviation inflation (not a correction
  to the vendor likelihood).
- `ape_onepass_likelihood_validation.py`: 128 correlated surrogate null
  and injected trials plus genuine vendor optimizer spot checks.
- `ape_observed_likelihood_profiles.py`: observed vendor likelihood profiles.
- `ape_analytic_kernel_crosscheck.py`: original pilot C comparison.
- `run_ape_axion_diagnostics.py`: one-entry-point 5/8 kHz workflow.

The 5 and 8 kHz original C pilot responses are 63.966461226 and
38.36500839 Hz^-1; the independent analytic responses were
63.9664611939 and 38.3650070861 Hz^-1.
For the 81-bin signal-weighted covariance diagnostics, the approximate
standard-deviation inflation ranges 1.78–1.82. This is **not** a
validated likelihood-ratio calibration or confidence-interval coverage.

## Reproduce on VIPER

From the repository working copy containing `APE_equivalent_external`:

```bash
cd /path/to/APE-analysis-code/APE_equivalent_external
export PYTHONPATH="$PWD"
ROOT=/ptmp/qarokn/APE_4p28h_PSD_OriginalVendor_PostLPSD_20261117
python -m unittest discover -s ape_post_lpsd/tests -p 'test_ape_diagnostic_integration.py' -v
python ape_post_lpsd/run_ape_axion_diagnostics.py \
  --root "$ROOT" \
  --output-dir "$ROOT/diagnostics/ape_repro_20261119_run01"
```

Use the Python environment with the original vendor dependencies installed.
**The output directory must not exist** because the workflow refuses to
overwrite data. It checks existing real X/Y spectra and fitted backgrounds.

Default **working** optical parameters used for reproducing the prior preview:
finesse_s=1500, tau_s=2.5e-9 s, tau1_s=1e-8 s,
rho_dm_gev4=3.1e-42 GeV^4, unit_conversion=1.519e24, and
g_ref=1e-10 GeV^-1. These are not a freshly certified APE calibration.
The reference signal PSD is `0.5*g_ref**2*sum_k(h_k*R_QWP(f_k)**2*K_D(j,k))`.
The original C output is the D PSD column (rad^2/Hz), not S.

The runner intentionally reproduces **only** the two benchmark masses,
81 bins around each, and separate X/Y diagnostic vendor likelihood calls.
It does not extrapolate or silently claim the full 5,035,956-bin signal
kernel. The window uses the tested *pilot* PSLL=238.13 dB; matching
production window provenance is still outstanding. The fixed laboratory
halo is a benchmark with no annual modulation, not a certified lab halo.

## Interpretation and next science gates

The injection validation used an intentionally strong signal (peak 2.5
times the fitted background); passing is an engineering smoke test, not
evidence of sensitivity at g_ref. The vendor likelihood remains an
independent-bin composite likelihood when evaluated on overlapping
Constant-Q bins. The 1.8x diagnostic inflation must not be applied blindly
to vendor sigma_mu. The observed profiles are uncalibrated.

To progress to a physical full-band scan: verify production window
parameters, calibrated complex two-QWP optical transfer and X/Y projection,
validated lab-frame halo, and bin/quadrature covariance and Monte Carlo
coverage. Do not interpret negative fitted signed mu as negative physical
power. Keep original vendor statistical functions unchanged.

## Original vendor fidelity

Every new module above lives under `ape_post_lpsd`. Existing
`scalar_dark_matter_lpsd/` and `LPSD/` are not modified by this change.
