# Analysis Contract: Unified APE Constant-Q Spectral Workflow

Date proposed: 2026-10-16
Status: PROPOSED - researcher approval of this specific contract required before scientific code implementation.
Scope: integration of existing, validated C-LPSD production method; not a new estimator or new analysis stage.

## REQUESTED
Integrate the previously completed APE X/Y Constant-Q production workflow into the existing `APE-analysis-code` feature branch as one clean generic program, runnable standalone under Linux (including Docker) or on VIPER via Slurm. Reuse existing approved APE ingestion, cleaning, calibration, calibrated HDF5 outputs, and the unchanged external Scalar-Dark-Matter-LPSD executable. Do not copy redundant historical scripts wholesale.

## CURRENT BEHAVIOR
- Branch `feature/ape-dxd-ingestion-v1` at `fd624298a9e4819fa4d27da06f7a590a6efb3a0d` has APE ingestion, cleaning, calibration, HDF5, tests and an older method-0 comparison pathway.
- Docker `lpsd-linux` computed X via shared FFT (C method 3) and sixteen method-2 Constant-Q batches and contains X outputs and execution scripts.
- VIPER computed Y with the same numerical settings via Slurm arrays, and combined its sixteen outputs using the external project's combination script.
- Mac contains the two completed 5,035,956-bin spectra, their PSD sum, and ASD comparison plot. These are read-only reference outputs.
- The existing `src/ape_c_lpsd/cli.py` deliberately refuses production real-data execution. This existing gate must not be silently bypassed.
- `configs/run.yaml` and `src/ape_c_lpsd/core.py` describe a distinct older method-0 estimator. Do not relabel it as Constant-Q.

## UNDERSTANDING
One spectral implementation handles either `rho_x` or `rho_y`. Local and Slurm backends differ **only** in process submission and resource handling; their C arguments and output-validation rules are identical. Production input is a genuinely continuous, sample-ordered, float64 HDF5 quadrature, not X/Y magnitude and not concatenated disjoint intervals.

## PLANNED IMPLEMENTATION
1. Add a project-owned recorded Constant-Q profile, with explicit frequency grid, C options, and exact sixteen batch index boundaries.
2. Add a small spectral adapter that generates the actual method-3 shared-FFT C invocation once per input quadrature and method-2 C invocations for each of sixteen existing frequency batches.
3. Add execution adapters: local sequential invocation and Slurm script/array generation. Neither backend changes numerical/scientific settings. No job submission or real-data processing until separately authorized.
4. Add the existing-output reader/combiner with strict count, finite/nonnegative PSD, monotonically increasing frequency, no missing/overlapping batch indices, and exact X/Y frequency-grid consistency checks. Retain C PSD values and the `N_avg` and kernel-bin fields.
5. Add synthetic/fixture tests for command parity, batch boundaries, C output format, safe paths, and combination. Compare against existing validated spectra only by read-only checks where accessible.
6. Provide one documented user-facing entry point and preserve the existing approval gate. Do not duplicate the Docker/VIPER channel-specific shell files in the final program.

## SCIENTIFIC DEFINITIONS AND UNITS
Historical completed profile:
- Input HDF5 dataset: `rho_x` or `rho_y`, processed independently; `f_s=20000 Hz`.
- Frozen input: 308,717,000 samples = 15,435.85 s; *historical* C command used `-e 15435.0`, under its inclusive-endpoint convention. Preserve exact historical settings and provenance; do not silently extend the analyzed duration.
- `f_min=65 Hz`, `f_max=10000 Hz`, `J=5,035,956`, approximate relative spacing `1e-6` (`Q≈1e6`).
- Shared FFT: `-h 3`; Constant-Q averaging: `-h 2`. No second incoherent power averaging.
- Common options: `-A 2 -b 0 -f 20000 -l 20 -r 0 -w -2 -p 238.13 -x 1 -E 10 -M 26 -R 1e-15 -T`. Batch `-N` is **global start index**, `-n` is the number of requested bins, `-J 5035956` is the full grid.
- Batch boundaries (17 endpoints for 16 batches):
  `[0, 3649978, 3996383, 4199059, 4342872, 4454427, 4545576, 4622643, 4689403, 4748290, 4800966, 4848618, 4892122, 4932141, 4969193, 5003688, 5035956]`
- Raw method-2 default format was `fDNbri` (six columns); historic combined outputs retain the first four (`f, PSD, N_avg, N_kernel`). The new parser must distinguish formats explicitly and must not mistake the older method-0 five-column `fDNRb` format for the method-2 output.
- If `rho` is a polarization angle measured in radians, one-sided PSD is `rad^2/Hz`, ASD `rad/sqrt(Hz)`; preserve existing unit conventions.
- X/Y PSDs may be summed only after independent estimation on an exactly matching grid. No time-domain magnitude `sqrt(rho_x^2+rho_y^2)`.

## RISKS AND CONTROLS
- **Scientific drift:** compare all generated C options and batch boundaries against the captured Docker and VIPER scripts.
- **Wrong method:** fail closed if method 0 is selected in the new Constant-Q profile.
- **Missing/duplicate bins:** check each expected batch index/count and combined 5,035,956-bin coverage.
- **Large-memory files:** stream output rather than loading millions of text rows and repeatedly stacking them.
- **Backend differences:** generate Slurm jobs with explicit dependencies and prevent accidental duplicate execution/overwrite.
- **Results preservation:** output to new run-scoped paths only; no changes to preexisting spectra, raw data, or cached FFTs.
- **External dependency:** require a user-provided compiled C executable; external C/Python source is read-only.
- **Scratch:** VIPER `/ptmp` is ephemeral and requires separate archival decisions; never treat it as the only permanent copy.
- **Real-data execution:** require explicit later permission after synthetic validation.

## WILL NOT DO
No edits to ingestion, cleaning, calibration, external Scalar-Dark-Matter-LPSD sources, historical method-0 behavior, or existing results. No new signal processing, filters, averaging, background fit, axion likelihood, candidate cuts, upper limits, package installs, container recreation, cluster submission, or GitHub `main` merge. Do not commit raw HDF5, FFT, batch, or spectrum data to GitHub.

## FILES THAT WOULD CHANGE
New files only (final filenames may be narrowed before coding):
- `configs/constantq_ape_reference.yaml`: frozen numeric profile
- `src/ape_c_lpsd/spectral/constantq.py`: commands and batch plan
- `src/ape_c_lpsd/spectral/combine.py`: strict existing-output combination
- `src/ape_c_lpsd/backends/local.py` and `slurm.py`: execution adapters
- `tools/run_constantq.py`: optional explicit user-facing interface, preserving real-data approval gate
- `tests/test_constantq_integration.py`: fixture/synthetic tests
- `docs/CONSTANTQ_INTEGRATION.md` and updates to scientific equation reference where needed
Existing `README.md` may receive a short, accurate entry. No other changes without separate approval.

## VALIDATION PLAN
First run tests against small synthetic outputs and inspect generated command lines without executing real-data jobs. Validate exact sixteen-bin partition and method flags, units/columns, no gaps/duplicates, exact channel-grid match, strict refuse-overwrite, plus read-only comparisons with previously validated results and existing tests. Report discrepancies; no scientific parameter changes to make tests pass.

## UNRESOLVED DECISIONS
- Exact user-facing command name and external executable discovery may be finalized without changing numerical defaults.
- Real-data run authorization, any alternative duration policy, and any additional low-frequency/statistical science remain out of scope.

## REAL-DATA GATE
All current APE real data and outputs remain read-only. Do not submit or run any C spectral jobs on real data under this contract without a separate explicit researcher go-ahead.

## STATUS
WAITING FOR RESEARCHER APPROVAL OF THIS WRITTEN CONTRACT BEFORE CODE IMPLEMENTATION.
