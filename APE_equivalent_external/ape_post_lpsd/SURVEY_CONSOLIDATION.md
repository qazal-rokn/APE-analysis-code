# Survey consolidation

The 64-mass pilot survey completed on VIPER with 59 complete fits, 5
partial-sigma fits and 0 failed mass hypotheses. These results are for
**selected frequencies, not a continuous upper-limit curve**.

Use `ape_summarize_band_survey.py` to combine read-only per-mass checkpoints
into one auditable CSV, JSON summary and PNG:
```bash
cd /ptmp/qarokn/APE-analysis-code/APE_equivalent_external
ROOT=/ptmp/qarokn/APE_4p28h_PSD_OriginalVendor_PostLPSD_20261117
PY=/ptmp/qarokn/ape-python-venv-20261023/bin/python
PYTHONPATH="$PWD" "$PY" -m unittest discover \
  -s ape_post_lpsd/tests -p 'test_ape_summarize_band_survey.py' -v
PYTHONPATH="$PWD" "$PY" -m ape_post_lpsd.ape_summarize_band_survey \
  --survey-dir "$ROOT/diagnostics/ape_band_survey_64" \
  --output-dir "$ROOT/diagnostics/ape_band_survey_64_summary_v1"
```
The output directory must not already exist. This tool **never reruns**
the C-LPSD, halo/QWP signal computation, background fitting or vendor optimizer.
It validates the manifest hash and the 64 X/Y checkpoints.

The CSV contains original vendor signed `mu_hat`, raw `lnL0/lnLmax`,
nominal vendor `sigma_mu` (blank where the original FinDiff failed),
and the nominal one-sided `q0`. They are *uncalibrated* estimates from
correlated native C bins. They are not significances, upper limits,
physical coupling bounds or a combined X/Y result.
