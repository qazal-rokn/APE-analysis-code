# APE candidate and C PS controller integration review, 23 October 2026

## Approved scope
Implement main controller interfaces, connect an externally verified PS-enabled C executable, and reuse external candidate-search methodology where scientifically portable. Do not generate simulated data, validate end-to-end inference, calculate upper limits, edit original vendor/C, or modify main.

## Original candidate code audit
Original source: AlexandreGoettel/Scalar-Dark-Matter-LPSD, `scalar_dark_matter_lpsd/lvk/get_candidates.py` and `stats.py`. Candidate q0 is computed as 2*(logLmax-logL0) if best-fit mu>0, zero otherwise. A Gaussian-tail significance threshold is corrected by a nominal number of searched frequencies. Above-threshold contiguous bins are clustered. Then the original source refits H1/L1 detectors individually and applies agreement/power tests. It uses LIGO H1/L1 labels, original scalar beta transfer functions, one-bin assumptions, and includes an interactive pdb breakpoint. The original cluster routine also omits a candidate cluster open at the final bin. Reusing the complete script unchanged is therefore neither operationally nor scientifically appropriate.

## APE adaptation
`analysis/ape_candidates.py` implements one-sided q0 candidate preselection with optional caller-provided effective-trial approximation. It clusters adjacent candidates, including the terminal cluster, and returns a provisional list. X/Y are NOT substituted for H1/L1; their noise correlation and signal response require explicit modeling and a later, separately approved consistency/veto scheme. Global p-values cannot be claimed without null-distribution calibration.

`analysis/ape_candidate_scan.py` evaluates a specified positive-coupling finite grid using the unmodified vendor `stats.log_likelihood`, with `mu=(g/g_ref)^2`, APE QWP optical response and normalized vendor local peak shape. Fitted background, skew residual model and C-bin windows must be supplied and validated. It is a diagnostic grid scan, not the vendor optimizer or a coverage-calibrated significance.

## C integration
`analysis/ape_pipeline_controller.py` builds the existing frozen historical C plan plus explicit opt-in `ps_optin.plan_ps_optin` using operator-provided executable and its exact SHA256. C output remains native C, and raw7 sidecars are separated without changing the original C estimator. `tools/plan_ape_ps.py` is PLAN-ONLY. Full C submission/execution/finalization orchestration is NOT implemented by this milestone; the external binary and its own SHA256 must exist at runtime.

## Validation and limitations
Added focused candidate/PS reader checks. No test run was executed against the final checked-out GitHub commit in this session. Statistical trial factors, PS leakage/window response, exact peak-shape provenance, RMS/X-Y physics, correlated C bins, and background-fit selection remain outstanding. Real data and published significance forbidden. README updated. External vendor, C sources/binaries and original measurements remain immutable.
