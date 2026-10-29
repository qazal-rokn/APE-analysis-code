# APE 40-minute five-signal simulation: automation contract

Approved by researcher in chat on 2026-10-30 for automatic stage-by-stage execution and export, subject to existing scientific validation gates.

## REQUESTED
Run the *same* 40-minute, 20 kHz, five-injection synthetic dataset through the established APE pipeline on the self-hosted Mac runner, save each stage's inputs, numerical results, settings, logs, hashes, plots and replot scripts, and return a clear status per stage. Minimize manual Terminal interaction.

## CURRENT BEHAVIOR
Stage 3 ten X/Y sparse native-C PS products passed. Stage 4 ten vendor spline/skew-normal PS background fits passed. Stage 1/2 reports are now available from the archived original-fixture source and record successful original-APE synthetic HDF5-adapter executions: calibration of 48,000,000 samples per channel with maximum 0 rad reference difference and cleaning of 48,000 blocks with zero rejected. These isolated results do not validate native DXD input, disturbed-data cleaning, or actual Stage-1-output-to-Stage-3 file transfer. Stage 5/6/7 have not executed on this fixture. Local 3200-Hz vendor fit chi2=1.6836004673503357 and all other local fits passed; see immutable local results for the numerical detail.

## UNDERSTANDING
Automation is not authorization to improvise the APE optical response, the signal power normalization, the native C PS units, frequency-bin likelihood dependence or observed upper-limit convention. Preserve the five-signal original fixture (fixed signal phases; deterministic expected power; randomized detector noise), X/Y separation and scientific provenance.

## PLANNED IMPLEMENTATION
Phase A: on a manually attributable push to this repository's default branch, run a read-only self-hosted-runner audit using fixed frozen hashes; archive *plot-ready* Stage 3/4 synthetic numerical results, standalone plots and plotting scripts as a GitHub Actions artifact; save an explicit full-stage status JSON and inventory. The full original HDF5 and FFT caches stay on the research Mac. No C/vendor execution in phase A.

Phase B: add an opt-in controller that verifies source/output SHA-256, automatically resumes from passed stages, executes only existing unchanged stage scripts under a new per-run folder, and stops on failures. First show validated output/frequency windows, then review Stage 5 choices below and approve their concrete Analysis Contract. Independent plotting must never rerun C or refit vendor statistics.

## SCIENTIFIC DEFINITIONS AND UNITS
Native C column 3 PS [rad²] for calibrated ellipticity; native PSD column 2 [rad²/Hz], not interchangeable. Original single-frequency-bin vendor likelihood differs from unvalidated correlated multi-bin halo-template likelihood. Keep q0/local-significance vs global trials separate and expected upper limits distinct from observed physical exclusions. Reference coupling and QWP transfer function must be declared and independently checked.

## RISKS AND CONTROLS
- Repo is private; self-hosted runner label is ape-simulation; use workflow_dispatch or tightly limited workflow push for read-only preflight, not pull_request execution.
- No untrusted workflow checkout or native executable substitution.
- Reports may label Stage 1/2 as *isolated synthetic HDF5-adapter passed* only when original result reports, all retained samples, full source checksums and the relevant saved outputs are verified; no claiming Stage 5+ as done before computation.
- No modification to the original C binary/source, vendor scientific algorithms, frozen APE calibration/cleaning, raw data or previous results.
- Refuse missing/mismatched HDF5 or binary hashes.
- Original data and large FFT caches must not be uploaded to GitHub. Staged synthetic plots/spectra can be captured in private-run artifacts only.

## WILL NOT DO
No real-DXD execution; no changes to science formulas; no pseudo-signal generation, PSD-to-PS conversion, fabricated completion markers, blind coherence/correlation corrections or physical exclusion claim.

## FILES THAT WOULD CHANGE
Dedicated audit script in `tools/`, read-only workflow in `.github/workflows/`, this contract. Later orchestration modules only after a reviewed Stage 5 implementation contract.

## VALIDATION PLAN
Run read-only Mac job; verify fixed hashes and all 10 native-C Stage 3 PS spectra; verify original Stage 4 pilot + nine backgrounds, fit quality, full CSV/NPZ and separate plots; save SHA-256 manifest. Status audit should show incomplete/unverified stages and exact scientific blockers. No changes to already accepted files.

## UNRESOLVED DECISIONS
1. Decide and document the single-bin original vendor profile likelihood versus scientifically validated multi-bin halo likelihood; selected frequency points and X/Y combination must be explicit.
2. Confirm the signal-power normalization and optics response for the five-signal fixture; synthetic injected nrad RMS does not automatically provide the physical g_{aγγ} reference response.
3. Validate independent observations / native C spectral correlations before assigning global significance or observed-limit coverage.
4. Choose and test negative MLE boundary treatment if observed limits are later needed.

## REAL-DATA GATE
No real APE physical inference authorized. All automated products must remain labelled synthetic, provisional, and local unless independently validated.

## STATUS
Phase A read-only GitHub automation authorized; Phase B scientific Stage 5–7 gated by documented, validated model choices. This contract does not certify the entire future analysis chain.

## STANDING USER SCOPE RULE: only this five-signal fixture

The researcher explicitly authorized proceeding independently for this **specific** existing five-signal synthetic fixture without repeated operational approvals. No unrelated GitHub files, branches, workflows, experiments, real data or scientific code are in scope. Only the simulation-specific workflow, automation/reporting wrappers, contract and saved synthetic results may be changed. Never replace, patch or rebuild the native C executable or original external vendor/APE scientific methods.

When any stage fails or an interpretation requires an undocumented assumption, examine the **original Göttel–Raymond 2026 paper** (DOI 10.1038/s41598-025-33428-2), the archived external vendor implementation and existing source manifests. Use the unchanged original method and preserve relevant source/code hashes. Do not invent alternative filters, reference couplings, signal templates, units, likelihood degrees of freedom or statistical coverage.

The approval covers operation and exact reuse of independently validated methods, **not** fabrication of missing scientific definitions. If a required parameter or scientifically valid method is unavailable in original sources or frozen contract, write a traceable `BLOCKED_NEEDS_SCIENTIFIC_DEFINITION` status and preserve all prior artifacts. Do not declare stages 5–7 passed by relabelling outputs. The user does not need to approve routine file archiving, reproducibility checks, job retries of unchanged methods or status reporting for this fixture.

A ChatGPT session cannot control the Mac or run continuously by itself; execution beyond the active chat requires the connected GitHub self-hosted runner to be online and an actual GitHub Actions workflow to invoke the declared unchanged stage executables. Do not claim that stages have been run without checking the returned Actions job logs and resulting artifacts.
