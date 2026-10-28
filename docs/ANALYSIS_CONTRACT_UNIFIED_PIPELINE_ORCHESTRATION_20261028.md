# Draft analysis contract: Unified APE pipeline orchestration (2026-10-28)

## REQUESTED
Reuse previously validated APE calibration, cleaning, native C-LPSD, original-vendor background fitting, likelihood, candidate and expected-upper-limit components. Connect tested stages into one automated, auditable pipeline; do not repeat past validation unnecessarily or modify authoritative scientific algorithms.

## CURRENT BEHAVIOR
tools/calibrate_frozen_segment_to_hdf5.py already uses the verified frozen-segment manifest and existing calibrate_chunk to write calibrated HDF5 rho_x and rho_y, plus i0_w and eta0. tools/run_ape_pipeline.py already executes the PS-enabled native C for X/Y, combines native PS sidecars, runs original vendor background fitting and a provisional finite coupling grid candidate scan. Its strict synthetic-only and halo provenance gates are intentional. Separately tested ape_vendor_fit.py and ape_vendor_candidate_report.py call the unchanged original vendor MLE, one-sided q0, candidate clustering and expected UL but are not called by the main runner. Bounded local native-C 5 kHz injection and matched-noise controls passed functional checks; physical amplitude recovery and coverage remain unresolved.

## UNDERSTANDING
Preserve all accepted prior validations within their tested scope. Integrate without modifying the original paper statistical assumptions, instrument calibration, C kernels or vendor likelihood. First integration must remain synthetic-only and single-bin, with provisional physical interpretation.

## PLANNED IMPLEMENTATION
1. Add explicit opt-in vendor-MLE mode to tools/run_ape_pipeline.py, retaining existing finite-grid diagnostic path.
2. After C PS and vendor background fits, prepare explicitly selected single-bin candidate windows with the existing vendor-shape APE interface; call existing fit_and_report_single_bin_candidates, which delegates unchanged original vendor MLE, q0, candidate grouping and expected-UL calculation. Reject unsupported multibin shapes in this mode.
3. Generate complete no-overwrite manifest and CSV: calibrated input provenance, C binary/version hashes and all effective flags, background fit quality, signed muhat, sigma_mu, q0, local significance, candidate clusters, expected UL and explicit failure statuses.
4. Test against existing immutable native-C outputs (no expensive native C rerun necessary), negative-input cases, and original vendor numeric parity. Optionally run one bounded native-C end-to-end synthetic confirmation.
5. Keep frozen-segment cleaning/calibration implementation unchanged; report its provenance and existing independent entry point but do not fake synthetic provenance for real data. Integrating automatic fresh .dxd ingestion/cleaning into a real-data workflow is a separate, separately authorized scientific milestone.

## SCIENTIFIC DEFINITIONS AND UNITS
rho_x/rho_y are calibrated radians. Native C direct PS is rad^2, not PSD rad^2/Hz. Vendor receives log observed and background PS, reference signal peak_norm, normalized one-bin shape [1], original skew-normal params. mu=(g/g_ref)^2. q0=2(Lmax-L0) if muhat>0, else zero. sqrt(q0) is local ranking statistic, not global significance. Vendor expected 95% UL is asymptotic sqrt(mu_UL), not observed upper limit or demonstrated 95% coverage. All physics configuration, g_ref and native-C settings must be explicit.

## RISKS AND CONTROLS
The multi-bin paper halo/C PS kernel, C response normalization, peak/RMS and X/Y physics, global trial factors, and APE 95% limit coverage are not yet scientifically closed. Vendor FinDiff sigma may fail: record invalid output explicitly, do not silently replace with a valid limit. C local grids require minimum 1100 background frequencies and physically valid sample windows. No arbitrary signal model or ad hoc normalization.

## WILL NOT DO
No modifications to main, original external vendor C/Python, frozen reference calibration or cleaning, raw .dxd, calibrated caches, past results, or synthetic provenance gates. No Python approximation of C-LPSD; no redefinition of likelihood or limits; no fake physical detection claims; no heavy Viper-sized run on Mac.

## FILES THAT WOULD CHANGE
tools/run_ape_pipeline.py; a new focused tests/test_ape_unified_pipeline.py; an explicit run config in configs/; README.md; docs/EQUATION_REFERENCE.tex; this contract. Further changes need specific prior approval.

## VALIDATION PLAN
Vendor optimizer and UL parity on existing test arrays; exact single-bin scan and terminal cluster behavior; reuse earlier C spectra as read-only test inputs; output manifest/provenance and failure gates; single bounded synthetic native-C confirmation only if needed. Report all mismatches between request, contract, implementation, validation and documentation.

## UNRESOLVED DECISIONS
Candidate scan index/range and blind-selection policy; how many vendor-sigma failures permit a run to complete; later real-data orchestration; physically correct C PS response and halo signal shape.

## REAL-DATA GATE
NOT AUTHORIZED. No passing real/hybrid observations through a synthetic-only gate by changing metadata.

## STATUS
DRAFT: AWAITING EXPLICIT RESEARCHER APPROVAL OF THIS ANALYSIS CONTRACT BEFORE SCIENTIFIC CODE CHANGES.
