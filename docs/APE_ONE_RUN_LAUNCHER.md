# APE one-command launcher: current verified scope

Entry point: `python -B -m tools.ape_one_run <mode> ...`.

The launcher delegates to the existing implementation; no additional FFT, interpolation, statistical estimator, dark-matter realization, or vendor function is invented. Only modes corresponding to implemented, distinct scientific contracts are available.

| Mode | Covered stages | Boundaries |
| --- | --- | --- |
| `synthetic-full` | Existing synthetic HDF5 -> verified native C power spectra -> original vendor background -> existing likelihood/candidate diagnostic | Requires existing explicit synthetic input/approved executable config; its original synthetic-only guards remain in force |
| `frozen-17freq` | Existing native C X/Y spectra -> deterministic synthetic halo power -> original vendor background -> likelihood -> local candidate clustering -> expected 95% limits -> regression reference | Reuses already computed C spectra, does not rerun C; frozen 12-bin/17-frequency benchmark only |
| `real-calibration` | Verified frozen real DXD -> sample-level RMS calibration -> HDF5 -> validation and provenance | Mac ARM64, dry-run unless `--execute-calibration`; no C or physical coupling inference |
| `real-calibration-reuse` | Non-mutating metadata provenance check of already completed calibrated HDF5 | No C or vendor analysis, optional costly complete SHA-256 recheck |

## Existing real data: reusable, do not recalibrate

```bash
cd "$HOME/Desktop/APE_Research/APE_Automatic_Orchestrator"
export PYTHONPATH="$PWD/src:$PWD"
export PYTHONDONTWRITEBYTECODE=1
PY="$HOME/Desktop/APE_Research/APE-data-analysis-c-lpsd/.venv/bin/python"
arch -arm64 "$PY" -B -m tools.ape_one_run real-calibration-reuse \
  --stage-result "$HOME/Desktop/APE_Research/APE_Automatic_Runs/ape-run-20261028-c0284-arm64-v2/stage_result.json"
```

This reuses the successfully calibrated 308717000-sample, 20000-Hz experimental segment with C_trans=0.0284 W/V, q=0.7 A/W, G=1e6 V/A and H2=0.7 (provisional measurement provenance). Its HDF5 has the recorded SHA-256 `1edd21da7af3cc3026c7eb3b5557b638a31f41e67fcef030328ecb922dadab42`.

## Synthetic 17-frequency regression mode

Supply the *existing frozen files* from the benchmark, plus a fresh non-overwriting output directory:

```bash
/work/ape-statistics-venv/bin/python -B -m tools.ape_one_run frozen-17freq \
  --rho-x /absolute/path/to/frozen_x.ps \
  --rho-y /absolute/path/to/frozen_y.ps \
  --model /absolute/path/to/frozen_deterministic_model.json \
  --config /absolute/path/to/frozen_postprocess_config.json \
  --output-dir /work/ape-postprocess-new-unique-name \
  --reference-report /absolute/path/to/frozen_17freq_summary.json
```

The `synthetic-full` mode takes the **existing synthetic-only JSON** used by `tools.run_ape_pipeline` and delegates without altering its inputs/approval checks.

## Real-data provisional full-chain entry point (29 October 2026)

`python -B -m tools.ape_one_run real-provisional-full --config configs/ape_real_one_command.example.json` now calls the actual existing frozen calibration and direct-C PS flow, original external vendor background fit/MLE/expected-limit functions, APE local candidate grouping and provisional final CSV/plot in a single invocation. Use `--dry-run` for a guarded preflight. The run configuration supplies frozen scientific settings and initial environment/input paths, **not manually connected intermediate files**. Outputs are non-overwriting; previously computed PS can be reused only with matching calibration SHA-256 provenance.

This is deliberately **not** a general raw arbitrary DXD analysis, a validated physical multibin axion-halo inference, observed exclusion limit, or global candidate significance. The current external vendor MLE bridge supports only the existing one-bin approximation; replacing that with a physical axion multibin response is outside orchestration scope and is not silently done. See `docs/APE_REAL_ONE_COMMAND_INTEGRATION_20261029.md`.

## Important outstanding integration

This is **one CLI launcher** for the already implemented workflows, **not yet a single uninterrupted real DXD -> native C -> measured-data candidates -> physical 95% limits pipeline**. A real-data-only C backend with a different non-synthetic authorization contract must be developed and validated; do **not** reuse `run_local(synthetic_approved=True)` to bypass this. Do not run the frozen synthetic 12-bin power model on experimental input and label the resulting coupling limits physical. Final QWP X/Y complex response, native C correlated residual coverage/global trials, and negative-MLE observed limits remain unresolved. The 256-trial coverage diagnostic applies only to its tested synthetic statistical model and is not an automatic per-run stage.

The existing original C and vendor programs, raw DXD, and prior benchmark/frozen reports must remain unchanged. Do not overwrite existing outputs. Both `main` and frozen reference products remain untouched.
