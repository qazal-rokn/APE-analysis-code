# Resumable provisional APE band-wide axion likelihood survey

`ape_band_likelihood_scan.py` extends the committed and previously tested
5/8 kHz halo + analytical Kaiser + QWP + **unchanged original vendor** likelihood
workflow to a reproducible logarithmically spaced selection of independent
axion rest-frequency hypotheses across the APE band.

**This is a selected-mass diagnostic scan, NOT a calculation at every one of
the 5,035,956 native C output bins and NOT a calibrated 95% physical upper limit.**
The original C and vendor statistics source files are unchanged.

## VIPER

From `/ptmp/qarokn/APE-analysis-code/APE_equivalent_external`:

```bash
export PYTHONPATH="$PWD"
/ptmp/qarokn/ape-python-venv-20261023/bin/python -m unittest discover \
  -s ape_post_lpsd/tests -p 'test_ape_band_likelihood_scan.py' -v
ROOT=/ptmp/qarokn/APE_4p28h_PSD_OriginalVendor_PostLPSD_20261117
PYTHONPATH="$PWD" /ptmp/qarokn/ape-python-venv-20261023/bin/python \
  -m ape_post_lpsd.ape_band_likelihood_scan \
  --root "$ROOT" --output-dir "$ROOT/diagnostics/ape_band_survey_64" \
  --count 64 --fmin 65.1 --fmax 8000 --start-index 0 --stop-index 4
```

Next chunks: set `--start-index 4 --stop-index 8`, etc, until 64.
Use a compute node or batch job for large surveys, not the shared login node.
`scan_manifest.json` holds the exact selected frequencies and config hash;
`mass_0000000/completed.json` is a successful checkpoint with X/Y fits.
Missing/failed masses have `error.txt` and can be safely retried.
Changing the configuration requires a different output directory.

Masses are log spaced, individual halo + C kernels occupy 81 native bins
around the rest frequency, and vendor background segments must cover the
whole 81-bin band. Boundary-straddling or bad-fit mass points are recorded as
failures, not silently passed or interpolated.

All reported likelihood statistics are **uncalibrated** owing to strongly
overlapping frequency bins (~0.84 lag-one correlation near 5/8 kHz).
The pilot Kaiser PSLL=238.13 dB, optical QWP amplitude/phase and X/Y
signal fractions, lab-frame halo and proper Monte Carlo coverage remain
provisional. A band survey does not resolve these scientific gates.

## Handling the 75.842844 Hz vendor FinDiff sigma failure

On the first VIPER chunk, 65.1, 70.266415 and 81.861824 Hz passed
X/Y. At 75.842844 Hz, the unchanged vendor optimizer returned finite
`lnL0`, `lnLmax`, and `mu_hat` for X and Y, but its
`stats.sigma_at_point` failed with `Spacing must be > 0` in Y.
The APE scan wrapper now records `mass_0000002/partial.json`
(`PROVISIONAL_PARTIAL_SIGMA_UNAVAILABLE`) and retains the existing
likelihood fit instead of treating this as a broken axion model.

`completed.json` is reserved strictly for points with both finite
vendor likelihood fits **and positive finite sigma_mu**. A partial
checkpoint is not counted as a complete likelihood result and is not
permitted to become a reported upper limit. Neither original vendor
`stats.py` nor `max_lkl_executable.py` is modified. Revisiting the
Fisher numerical issue and correlated-bin statistical coverage is a
separate validated work item. Failed mass-point input/output files
remain read-only; the new wrapper consumes the original `likelihood.json`
on resume without rerunning the C estimator, halo/QWP preparation, or
vendor optimizer.

To reclassify the existing partial point and continue with new masses:

```bash
cd /ptmp/qarokn/APE-analysis-code
git pull --ff-only
cd APE_equivalent_external
ROOT=/ptmp/qarokn/APE_4p28h_PSD_OriginalVendor_PostLPSD_20261117
PYTHONPATH="$PWD" /ptmp/qarokn/ape-python-venv-20261023/bin/python \
  -m unittest discover -s ape_post_lpsd/tests \
  -p 'test_ape_band_likelihood_scan.py' -v
PYTHONPATH="$PWD" /ptmp/qarokn/ape-python-venv-20261023/bin/python \
  -m ape_post_lpsd.ape_band_likelihood_scan \
  --root "$ROOT" --output-dir "$ROOT/diagnostics/ape_band_survey_64" \
  --count 64 --fmin 65.1 --fmax 8000 \
  --start-index 2 --stop-index 8
```

## Cross-segment backgrounds (88.358478 Hz)

The pilot scan detected that the original 81-bin native-frequency template
may straddle a fitted-background segment boundary. Previously the APE wrapper
required one segment to cover the whole frequency window and aborted.

The APE adapter now **retains all bins** and assigns each its correct original
vendor-fitted spline and skew-normal parameters. In a single-segment window
it preserves the original one-row layout; at a boundary it reshapes all bins
into 81 one-bin rows, using the same signal amplitude `mu` for every row.
This calls the original vendor `calc_max_lkl` and `stats.log_likelihood`
unchanged and preserves the halo PSD, spectrum normalization and native grid.
It rejects uncovered bins, ambiguous interior overlap, and chi-square failures.
The original vendor source is not modified. Like all other pilot results, the
overlapping-bin composite likelihood is NOT statistically calibrated.

The 88.358478 Hz failure has existing `halo.npz` and `qwp.npz` inputs,
which the resumable scan will reuse. Run the integrated tests before retrying:

```bash
cd /ptmp/qarokn/APE-analysis-code
git pull --ff-only
cd APE_equivalent_external
PY=/ptmp/qarokn/ape-python-venv-20261023/bin/python
ROOT=/ptmp/qarokn/APE_4p28h_PSD_OriginalVendor_PostLPSD_20261117
PYTHONPATH="$PWD" "$PY" -m unittest discover -s ape_post_lpsd/tests \
  -p 'test_ape_segmented_vendor_background.py' -v
PYTHONPATH="$PWD" "$PY" -m ape_post_lpsd.ape_band_likelihood_scan \
  --root "$ROOT" --output-dir "$ROOT/diagnostics/ape_band_survey_64" \
  --count 64 --fmin 65.1 --fmax 8000 --start-index 4 --stop-index 12
```

This preserves the scan manifest/hash and avoids recomputing completed points.
