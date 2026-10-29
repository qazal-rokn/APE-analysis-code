# APE verification baseline freeze — 2026-10-29

**Source:** `development/likelihood`. **Scope:** preservation of existing verified code and evidence; does not assert full end-to-end scientific validation, GitHub CI pass, or real APE discovery/exclusion. **Frozen by:** GitHub snapshot branch `freeze/ape-verified-baseline-20261029` created after this documentation commit. Do not commit additional work to the freeze branch. Future development remains on `development/likelihood` or a new feature branch.

## Preserved stage and evidence inventory

| Stage | Verification evidence and scientific scope | Baseline status |
| --- | --- | --- |
| DXD ingestion / selection and cavity/eta0 cleaning | Prior approved end-to-end ingestion and frozen cleaning tests; `docs/ANALYSIS_CONTRACT_APE_CLEANING.md`, `tests/test_cleaning.py`, `tests/test_ingestion_segments.py`; earlier longest-segment verified on actual acquisition | **Previously verified**; do not change implementation/rules |
| RMS calibration of separate `rho_x`, `rho_y` | `docs/ANALYSIS_CONTRACT_APE_CALIBRATION.md`, `tests/test_calibration.py`; prior forward synthetic recovery; later 0.0284 W/V real run separately provenance-bearing | **Previously verified as calculation/software**; H2 and physical optical response remain provisional |
| Original native C-LPSD, power units and validated frequency grid | `docs/VALIDATION_MILESTONE_LPSD_2026-10-28.md`, `docs/APE_SPECTRAL_UNITS_20261029.md`, C synthetic PS recovery tests. Real X/Y native PS each 5,035,956 rows from 65–10,000 Hz exist on Mac | **Previously verified as estimator/PS processing**; original vendor/C unchanged |
| Vendor background fitting | Original functions executed on synthetic pilot (4/4), 314 saved X real-data segments with zero recorded failures, `tests/test_ape_vendor_fit.py`; residual diagnostics showed deviations from selected model | **Software stage executed/validated in scope**, NOT proof of full-band noise distribution |
| Deterministic axion expected-power template | `docs/APE_DETERMINISTIC_SIGNAL_SIMULATION_CONTRACT_20261029.md`, `docs/APE_MANDATORY_DETERMINISTIC_DM_SIGNAL_SIMULATION_RULE_20261029.md`, 5-kHz frozen reference regeneration | **Passed deterministic paper-assumption / numerical checks**; QWP response and X/Y physical partition not finalized |
| Vendor likelihood/MLE and local candidates | Mac paper-consistent 17-frequency test: 34/34 successful after adapter input representation fix; injection `mu=4` recovered `mu_hat=3.5767197885` at 5000 Hz, 8.1255 local sqrt(q0); one injected cluster, zero control clusters; `tests/test_ape_vendor_candidate_report.py` | **Synthetic single/limited-frequency software checks passed**; local only |
| Vendor asymptotic *median expected* 95% g limit | Mac 34/34 formula-parity results; max relative discrepancy `4.498908935979278e-12`; frozen 17-frequency benchmark expected `g_UL(5000 Hz)=7.298703022667408e-11 GeV^-1` applies to its separate reference fixture | **Formula and numerical consistency passed**; not observed/coverage guarantee |
| Original-paper-style one-bin profile likelihood | Mac 1000 hypotheses × 512 realizations each, paper-style skew-normal synthetic log-power residuals: original vendor parity 80 checks `5.68e-14`; IID null 5.28%, across-frequency correlated null 5.25% vs nominal 5%; deterministic `mu=4` recovery 4.0171/4.0191, q_mu coverages 95.19%/95.25% | **Controlled one-bin local model passed**; one-bin conclusion must not be applied to correlated *multi-bin joint* inference |

## Diagnostics NOT validated or accepted as final science

- Real native-X 65–89 Hz bin correlation (lag1 approximately 0.795), and diagnostic multi-bin score inflation (~3.35) highlight possible problems **if** several correlated bins are multiplied into the same joint likelihood. This is NOT evidence that real APE peaks are false, and does NOT invalidate a correctly specified one-bin local likelihood.
- Multibin 2026 halo-template local significance and its dependence/coverage: not approved. Simulated deterministic *power* injection is after C and does not validate a time-domain axion propagation through C.
- True global candidate significance and trials correction: not approved.
- Full-band real APE noise model, X/Y independence, and publication-grade 95% **observed** upper-limit coverage (including signed/negative MLE and boundaries): not approved.
- Physical QWP transfer function, X/Y signal-power fraction and measured `H2` transfer: provisional.
- The one-command real-data orchestration under `tools/ape_real_one_command.py` is integration code with a provisional single-bin bridge; NOT independently demonstrated as a fully successful user-Mac real-data one-command output. No claim of completed end-to-end production validation.

## Permanent invariants

- No modification to original C-LPSD or vendor source/binaries, original cleaning/calibration calculations, real raw/cached data or prior benchmarks.
- `C_trans=0.0284 W/V`, `q=0.7 A/W`, `G=1e6 V/A` and sample-level `eta0(t)`, separate RMS X/Y; record actual run configuration and `H2` status.
- Native PS is `rad^2`; do not substitute PSD, invent coherence multipliers, or apply `N^{-1/4}` as an after-the-fact correction.
- Before *every* new signal simulation, review the original 2026 Göttel–Raymond paper + vendor implementation. In the paper approximation expected DM power is deterministic: no random DM phase/amplitude. Detector noise may be randomized.
- Distinguish local one-bin method from an unvalidated multi-bin likelihood; distinguish expected from observed limits; never claim global significance without its trials calibration.
- Each future validation requires separate test artifacts, immutable inputs, tracked revision/config hashes, numerical assertions, and a documented acceptance decision.

## Evidence provenance limits

The quantitative Mac runs above are based on user-provided Terminal output from this project conversation, not uploaded GitHub CI artifacts. The frozen documentation is an **audit record** and branch snapshot, not independent re-execution of all local tests or a certification of the complete pipeline.

## Exact next step

On an **active development branch**, perform one end-to-end *single-bin pilot* from existing, checksum-verified `rho_x.npy` / `rho_y.npy` native C spectra to cached vendor background, original one-bin local likelihood/MLE, candidate clustering, and median expected limit at a **small selected frequency set**. Reuse the existing native PS and background files; do not recompute C or calibrate raw DXD. Verify provenance, file shape, 65–10000 Hz grid, native PS units, source hashes, and one-command output parity against individually checked functions. No real detection/exclusion claims. Before production, separately validate the physical coupling model and actual observed-limit/global-significance conventions.
