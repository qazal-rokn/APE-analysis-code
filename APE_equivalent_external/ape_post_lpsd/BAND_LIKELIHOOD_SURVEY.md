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
