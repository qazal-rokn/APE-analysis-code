# Analysis Contract: Opt-in multi-bin APE forward template (DRAFT, 2026-10-28)

## REQUESTED
Implement a physically justified multi-bin APE axion signal-power model in the existing unified synthetic-only C-LPSD -> original vendor likelihood pipeline, avoiding repetitive C normalization tests or extra ad hoc calibration factors.

## CURRENT BEHAVIOR
- `src/physics/qwp_signal.py:qwp_signal_power_template` already computes a provisional `T[q,j] = g_ref**2/2 * alpha[q] * sum_k h[k] R_QWP(f[k])**2 W[q,j,k]` using caller-supplied halo weights, quadrature power fractions and C window power responses. Inputs are PS (rad^2).
- `src/ape_c_lpsd/analysis/ape_likelihood.py:prepare_ape_likelihood` already factors this array for the unchanged original vendor `stats.log_likelihood` and logs observed/background PS.
- `src/ape_c_lpsd/analysis/vendor_shape_ape.py` accepts an external vendor shape already including C window effects, but does not generate such effects.
- `tools/run_ape_pipeline.py` supports a `finite_grid` mode and `vendor_mle_single_bin` mode. The latter, `ape_vendor_fit.py`, and `ape_vendor_candidate_report.py` explicitly reject >1-bin MLE. No full multi-bin physical MLE/UL integration has been authorized.
- Completed coherent tone C tests at Q=1e6 show centered K=1, half-bin K≈0.9235237, Kaiser analytic parity <2e-6, and validated 7x7 overlapping-bin response. These do not validate the physical halo or joint-bin likelihood.

## UNDERSTANDING
The model must represent expected **native C-estimator bin powers** after window response, not physical input halo weights alone; **no extra K correction** on observed PS, coupling, or native C. The 2026 Göttel-Raymond paper concerns scalar DM in LIGO, not the APE pseudoscalar optical response; its halo and covariance assumptions must not be transplanted as established APE physics.

## PLANNED IMPLEMENTATION (requires researcher approval)
1. Create an explicit, opt-in `physical_multibin_template` synthetic-only mode. Reuse existing original compiled C X/Y, unchanged original vendor background and unchanged original `stats.log_likelihood`.
2. Load **externally supplied and documented** intrinsic-field frequency support and nonnegative normalized weights (halo source and convention mandatory), and a precomputed C-window power matrix `W[q,j,k]` with unambiguous frequency grids/units; refuse mismatches. No default astrophysical halo, guessed per-bin signal powers, or silent re-normalization/interpolation.
3. Generate PS `T[q,j]` only with the current `qwp_signal_power_template` and separately declared QWP inputs, `g_ref`, and X/Y fractions. Record full source paths, hashes, effective settings, and template PS; use `prepare_ape_likelihood`/original vendor log likelihood. Preserve the existing single-bin and finite-grid modes without behavior change.
4. Provide a scientifically labeled **multi-bin likelihood evaluation** for explicitly chosen candidate windows and user-provided coupling values. Retain all candidate results as *diagnostic*; do **not** silently reuse the single-bin vendor MLE or expected UL. Extending the vendor MLE to multi-bin with valid uncertainty and expected UL is an additional decision gated by review of source optimizer assumptions and C-bin covariance.
5. Add focused validation: exact input/shape/unit/no-overwrite guards, measured centered and off-grid window parity with archived C outputs, mathematical reconstruction `peak_norm * peak_shape == T`, and unchanged original-vendor likelihood parity. No repeated heavy C execution.
6. Update the analysis contract and equation reference; clearly mark inference **provisional, synthetic-only**, without global significance or confidence-coverage claims.

## SCIENTIFIC DEFINITIONS AND UNITS
`rho_x/rho_y [rad]` -> native C `P[q,j] [rad^2]`.
`T[q,j] = (g_ref^2/2) alpha_q sum_k h_k R_QWP(f_k)^2 W[q,j,k]` in rad², with `sum_k h_k=1`, **no forced sums of W**. The source of X/Y fractions, optical QWP amplitude and complex phase, halo weights, C window handling and inter-bin correlations must be documented; mathematical implementation alone is not physical approval. Vendor model remains `B[q,j] + mu*T[q,j]`, `mu=(g/g_ref)^2`.

## RISKS AND CONTROLS
- Overlapping C bins are correlated; treating their likelihood product as fully independent requires explicit justification or a separately calibrated approximation. No global significance or 95% coverage claims.
- The supplied experimental 7x7 matrix is a narrow 5 kHz coherent-tone benchmark; do not silently extrapolate to other frequencies, window settings, or axion lines.
- Original scalar-LIGO peak-shape and transfer function must not be claimed to be an independently validated pseudoscalar APE model.
- Input provenance, physical phase, quadrature response fractions and coupling are still provisional.

## WILL NOT DO
No modifications to external C, original vendor, raw/frozen data, existing APE calibration, `main`, original outputs, run modes or original statistical assumptions; no independent random per-bin signal-power draws, no ad hoc extra normalization, no automatic real-data run; no broad simulation campaign.

## FILES THAT WOULD CHANGE
`tools/run_ape_pipeline.py` (opt-in mode only), a new opt-in data-loading/template adapter under `src/ape_c_lpsd/analysis/` if necessary, focused new tests under `tests/`, `docs/EQUATION_REFERENCE.tex`, and README or a focused usage document. `src/physics/qwp_signal.py` and original vendor are reused unchanged unless independently justified and separately approved.

## VALIDATION PLAN
Static contract/shape/unit tests and synthetic parity with archived Q=1e6 C and analytic Kaiser results. One acceptance command for focused tests. No new C-LPSD normalization runs.

## UNRESOLVED DECISIONS
- Which physically sourced intrinsic APE axion frequency/halo weights to supply (2026 scalar-paper assumptions are not automatically the physical APE halo).
- Whether an opt-in **diagnostic multi-bin likelihood evaluation** is acceptable first, while fully statistical multi-bin MLE/UL waits for covariance/model approval.
- Provenance for X/Y response fractions and QWP transfer parameters.

## REAL-DATA GATE
NOT AUTHORIZED; synthetic-only input checks remain mandatory.

## STATUS
**DRAFT, awaiting explicit researcher approval. Scientific code implementation not started.**
