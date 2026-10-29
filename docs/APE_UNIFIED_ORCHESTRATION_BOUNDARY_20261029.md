# APE unified orchestration boundary (2026-10-29)

`python -m tools.ape_unified_orchestrator --config <run.json>` invokes, in fixed order, eight **unmodified** configured stage executables: ingest/clean/select, calibrate, native C PS, background fit, vendor likelihood/MLE, candidates, expected upper limits, and final reporting.

Modes are `reuse` (pre-existing SHA-256-verified file) and `run` (exact existing executable argv and fresh output file). No shell commands, synthetic provenance bypass, implicit PSD-to-PS conversion, coherence multiplier, or changes to scientific formulas. The program verifies outputs, fails closed, and writes a nonoverwriting audit record. `--preflight` validates before executing.

## Important scope

This is a **generic, fail-closed orchestration layer**, not a completed zero-configuration real-DXD-to-physical-limits program. It requires a run JSON with all stage commands and artifact paths. Current `ape_one_run` implements distinct workflows; it does not yet supply scientifically authorized real-data C+vendor multibin or observed-limit executables. Do not fabricate these by relabeling synthetic inputs. The only allowed final label is provisional until the physical normalization, dependence/coverage, and observed-limit conventions are verified.

Run definition (each of eight stages):

```json
{"mode":"reuse","outputs":[{"path":"/absolute/existing.file","sha256":"64-hex-digest"}]}
```

or

```json
{"mode":"run","argv":["/absolute/python","/absolute/unchanged_stage.py","--existing-option","value"],"outputs":[{"path":"/absolute/new-output.file","sha256":null}]}
```

Run config must contain `schema_version: 1`, `run_dir` and a `stages` mapping with all eight exact stage names. Run directories and output artifacts must be new for each execution. The original external vendor/C files, frozen calibration and real DXD inputs must not be overwritten.
