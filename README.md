# APE Analysis / C-LPSD Project

This repository is evolving from the original compiled-C LPSD comparison into
the auditable APE analysis pipeline. Scientific changes remain approval-gated:
inspect first, agree the contract, implement only the agreed scope, validate,
then report.

## Canonical stage 1: Dewesoft .dxd ingestion

The first canonical pipeline block is now defined as:

```text
Dewesoft .dxd
    |
    v
explicit five-channel mapping
    |
    v
APEDXDData
```

It reads exactly:

- first-harmonic X quadrature;
- first-harmonic Y quadrature;
- second-harmonic X quadrature;
- second-harmonic Y quadrature;
- cavity transmission.

No calibration, cleaning, filtering, interpolation, resampling, trimming, gap
filling, segment selection, spectral estimation, or X/Y combination is allowed
inside this stage.

Copy `configs/dxd_ingestion.example.yaml` to a run-specific YAML file and set
the actual .dxd path, Dewesoft channel names, and run-specific calibration
parameters. Then inspect the file with:

```bash
python tools/ingest_dxd.py --config /path/to/run_config.yaml
```

The wrapper prints an ingestion summary only. It does not write calibrated data.

The run YAML already records the later calibration inputs
(C_{\mathrm{trans}},q,G,H_2), but stage 1 does not apply them. In particular,
`second_harmonic_pd_response` is (H_2), the measured attenuation/response of
the photodiode/electronics chain at (2f).

See `docs/DATA_FORMAT.md` for the raw-data contract.

## Existing C-LPSD comparison work

The repository also retains the earlier approval-gated infrastructure for
comparing the established APE Python/SpiCyPy LPSD workflow with the compiled C
`lpsd-exec` backend.

Controlled synthetic cache/transport validation exists; the older real-data
C-LPSD command remains deliberately refused until its separate scientific
contract is approved.

Reference inputs are read-only:

- `/Users/user/Desktop/Code/Projects/APE-data-analysis`
- `/Users/user/APE-main-code-inbox`
- `/Users/user/Downloads/Scalar-Dark-Matter-LPSD-1.0.0/LPSD`
- `/Users/user/APE_Local_Data/After_wedged_EOM/AXION_COHERENCE_CLEANED_LPSD/CALIBRATED_CACHE`

The project must never write into those locations.

## Layout

```text
configs/           run/backend settings
docs/              contracts, method notes, run records, changelog
notebooks/         optional presentation notebooks; no authoritative logic
results/           generated artifacts; source data never stored here
src/               ingestion, calibration-review, spectral/backend modules
tests/             synthetic and format-validation tests
tools/             thin, auditable command-line entry points
```

## Planned downstream pipeline

The intended long-term flow is:

```text
.dxd reader
    -> calibration
    -> LPSD
    -> background model
    -> APE axion signal model
    -> likelihood
    -> candidates / upper limits
```

Those later stages are separate scientific changes and are not implied by the
stage-1 reader.

## Existing ingestion / segment-review infrastructure

The older `src/ape_c_lpsd/ingestion/` chronology and quality-review modules
remain available as a reference/cross-check path. They reuse frozen reference
APE quality and calibration functions, split at invalid or unverified
boundaries, and can freeze the longest validated segment. They are not silently
invoked by the new canonical .dxd reader.

See `docs/INGESTION_AND_SEGMENT_SELECTION.md` and
`docs/INGESTION_ARCHITECTURE.md` for that legacy/reference workflow.

## Existing calibrated quadrature review

The review-only `calibrated_review` stage applies the frozen reference scalar
calibration denominator independently to the two recorded voltage quadratures.
It does not combine X/Y and remains a reference implementation until the new
APE calibration contract is frozen.

## Current likelihood development (23 October 2026)

**Repository policy:** `main` is the established APE baseline. The current
unfinished likelihood work is on `development/likelihood`; no new branches.
Update this README and the equation reference whenever pipeline code changes.

### Implemented on development/likelihood

- The APE-owned adapter `src/ape_c_lpsd/analysis/vendor_statistics.py` now
  preserves the original six-argument vendor likelihood contract:
  log observed/background (M,N), nonuniform `peak_norm` broadcastable to
  (M,N), normalized `peak_shape` (N,), and skew parameters (M,3).
  `factor_power_template` converts a reference-coupling PS template into
  the original norm/shape convention without changing vendor statistics.
- `src/physics/qwp_signal.py` provides the **provisional** QWP Case-2
  finesse-dependent optical response and a windowed PS template.
  Finesse, tau, tau1, dark-matter density, unit conversion, reference
  coupling, halo weights, C-LPSD window response and quadrature fractions
  are **mandatory caller inputs**, never embedded physical defaults.

### Not yet validated or integrated

The QWP amplitude/peak-to-RMS conversion and actual X/Y signal responses
are unvalidated. The halo line shape and actual C-LPSD PS kernel must be
provided by an independently verified model, not invented here.
The optical response has **no additional sqrt(f/Q_a)** integration factor.
Native PS has units rad², distinct from PSD in rad²/Hz.
The opt-in C direct-PS extension is still separate from this repository.
This is **not yet a completed or scientifically validated real-data
axion coupling/upper-limit pipeline**. Original external vendor, compiled
C, raw data and calibration remain unchanged.

Minimal interface tests live in `tests/test_vendor_statistics_adapter.py`
and `tests/test_qwp_likelihood_contract.py`.
See `docs/ANALYSIS_CONTRACT_QWP_LIKELIHOOD_20261023.md`
and `docs/EQUATION_REFERENCE.tex` for assumptions.

## Synthetic end-to-end integration: milestone 1 (23 October 2026)

**PS transport integration committed, execution not yet performed.**
The existing six-column C command planner is preserved. New
`src/ape_c_lpsd/spectral/ps_optin.py` provides an opt-in deep-copy wrapper
that requires an actual executable path and an independently supplied exact
SHA-256 of the experimental **PS-enabled** C binary. It changes only the
executable, C format (`fDSNbri`) and batch output suffix (`.raw7`).
The original C numerical flags and frequency partitions are unchanged.

The validated streaming converter reads the native C seven-column output
`frequency, PSD, PS, Navg, Nkernel, real, imaginary`. It preserves the
original numeric tokens in separate legacy six-column `fDNbri` products
and two-column `frequency,PS` sidecars. It enforces row counts, increasing
frequencies, finite values, nonnegative spectra, and no overwrite.
No spectral-estimation computation was replaced by NumPy, FFT or Welch.

**Status boundary:** The VIPER PS-enabled C executable and its complete
archived source are not available inside this session; the module has NOT
executed that binary. The archive SHA-256 is **not** the executable SHA-256.
This stage adds the opt-in transport interface only. End-to-end
synthetic-to-physical-likelihood testing and actual C PS execution remain
pending. The external C/vendor code and raw/calibrated data are untouched.
See `docs/PS_OPTIN_SYNTHETIC_INTEGRATION.md`.

Minimum local adapter check (no external C binary required):
```bash
PYTHONPATH=src python -m unittest tests.test_ps_optin
```

## Draft integration: C interface, halo and likelihood (23 October 2026)

On **development/likelihood** only:
- `src/ape_c_lpsd/spectral/ps_optin.py` binds the C plan to a caller-provided
  verified PS-enabled executable path and exact executable SHA-256; the
  compiled binary itself is **not** embedded or connected to a running
  machine in this repository. Original compiled C and original plan untouched.
- `src/physics/axion_halo.py` provides an optional **unboosted Maxwellian
  benchmark**, with explicitly user-provided rest frequency, speed scale,
  light speed, frequency bin edges and minimum captured probability.
  This omits Earth's/laboratory motion, halo escape cutoff and modulation;
  it is not a final halo model for physical axion limits.
- `src/ape_c_lpsd/analysis/ape_likelihood.py` creates a PS [rad²]
  physical QWP signal template from caller-provided halo weights and actual
  C-LPSD power-kernel response, factors it for the unmodified vendor
  likelihood, and maps coupling g [GeV^-1] to
  `mu=(g/g_ref)^2`. The vendor fitting and statistical functions stay external
  and unchanged. No new integration-time factor is inserted.
- `tests/test_ape_likelihood_wiring.py` covers the required minimal
  interface and normalization contracts. No simulation generation or
  end-to-end execution was performed in this milestone.

**Blocked for a genuine run:** the actual validated PS C executable on
Mac/VIPER and its own SHA-256; the actual C-LPSD PS window kernel;
experimentally justified X/Y response and peak/RMS convention; a
lab-frame halo model and suitable independent noise statistics.
No likelihood confidence limits or scientific detection claims are
authorized by this wiring alone. All calibration constants are caller inputs.

To run the focused interface tests in a checked-out environment:
```bash
PYTHONPATH=src python -m unittest discover -s tests -p 'test_ape_likelihood_wiring.py'
```

Continue work on `development/likelihood`. Keep `main` unchanged
until acceptance; document every pipeline modification here.

## Exact external peak-shape interface (23 October 2026)

`src/ape_c_lpsd/analysis/vendor_shape_ape.py` now directly accepts the
**original external scalar-field program's normalized peak-shape .npy**.
The input is copied unchanged: no replacement Maxwellian model, no
interpolation, and no silent renormalization. It rejects non-unit-sum,
nonfinite, negative or different-length shape arrays. This is an
alternative to the experimental benchmark in `src/physics/axion_halo.py`,
not an additional shape factor.

For each detector observation `q` and native Constant-Q frequency `f_j`,
the vendor receives:
- `peak_shape[j] = h[j]`, unchanged and summing to one;
- `peak_norm[q,j] = (g_ref * R_APE(f_j))**2 * alpha[q] / 2`;
- `mu = (g_agg / g_ref)**2`;
- log observed and background **PS**, with units rad².

This preserves the vendor's original
`B[q,j] + mu * peak_norm[q,j] * peak_shape[j]` structure; the
LIGO-specific scalar beta is **not** reused. The original
`stats.log_likelihood` is unchanged.

**Scientific guard:** This direct-shape interface is only correct if
the supplied peak-shape weights already describe the actual native
C-LPSD *observed* frequency bins after their window response. The
external shape file's length and normalization alone cannot prove this;
the original file's grid definition and kernel handling must be checked.
If the shape is defined on intrinsic-field bins instead, use the
general `ape_likelihood.prepare_ape_likelihood` interface with an
independently validated C-LPSD kernel and halo weights.
Peak/RMS and X/Y fractions remain injection-validation gates.

Focused unit-test file:
`tests/test_vendor_shape_ape.py`. No end-to-end run, original vendor
execution, or C-binary execution is claimed here. The full validated
PS-enabled C binary remains operator-supplied at runtime, guarded by
the opt-in executable hash.

## External scalar-model audit: sliding peak shape (23 October 2026)

The original external `max_lkl_executable.py` establishes that
`peak_shape` is a **normalized template applied to each consecutive
candidate window**, not an array covering the complete C-LPSD frequency
grid. For candidate starting index `i`, it selects the consecutive
`N=len(peak_shape)` bins `i:i+N` and passes the same shape to
`stats.log_likelihood`. The project does not contain a checked-in
`.npy` peak-shape example or its generation code. When the external
shape is omitted, the vendor uses the one-bin default `[1]`.

The APE adapter now accepts a normalized shape independently of the
full spectrum size and requires matching frequency-bin count only
**when evaluating the selected candidate window**. No LIGO-specific
frequency grid is assumed. We keep vendor shape values unmodified and
replace its scalar/LIGO amplitude response by APE's QWP response.

**Open provenance requirement:** to reproduce their *multi-bin*
physical signal shape precisely, obtain the actual user-supplied .npy
and its creation method, if one was used. Matching the same frequency
indices/Constant-Q settings and estimator power response remains to
be verified on synthetic data; this is not permission to modify the
signal shape speculatively.

## Candidate pipeline integration (23 October 2026; no upper limits)

**External source audit:** `scalar_dark_matter_lpsd/lvk/get_candidates.py`
uses a one-sided `q0 = 2(logL_best-logL_null)` for positive `mu`,
a trial-factor threshold, clustering, and H1/L1 detector consistency.
**The script is not portable unchanged:** it uses H1/L1 interferometer
identifiers, scalar-specific `get_beta`, a fixed one-bin candidate
assumption, hard-coded optimizer scan settings, and a live
`pdb.set_trace()` breakpoint. Its original terminal-cluster routine also
does not append a still-active final cluster. X/Y in APE are quadratures
of one experiment, not independent LIGO detectors.

**Implemented (development branch):**
- `analysis/ape_candidates.py`: vendor-inspired one-sided q0,
  optional *assumed* effective-trial threshold and contiguous candidate
  clustering; no H1/L1 consistency veto and no physical detection claim.
- `analysis/ape_candidate_scan.py`: finite positive-coupling grid scan
  using the unchanged external vendor likelihood and APE detector response.
  Produces best-grid coupling, null/best log likelihoods and candidate list.
  This is diagnostic ranking, not a calibrated significance or upper limit.
- `analysis/ape_pipeline_controller.py`: guarded native C PS command-plan
  constructor, 7-column batch sidecar finalization and two-channel
  candidate-window PS reader with exact frequency alignment.
- `tools/plan_ape_ps.py`: plan-only CLI requiring the PS-enabled C
  executable and its own SHA-256. It does **not** run C or submit jobs.
- `tests/test_ape_candidate_pipeline.py`: focused candidate clustering
  and PS alignment checks.

**Still missing for end-to-end execution:** the operator's actual verified
PS-enabled C executable; the original vendor peak-shape .npy file if
multi-bin signals are intended; a declared procedure for the native
C frequency-window response; the real-vendor runtime; a fully automated
16-batch execution/finalization pipeline and complete run manifest,
candidate CSV/plots, end-to-end synthetic tests and calibrated trials.
The existing source has a protected route for these operations, but a
single input-to-candidate execution command is **not yet complete**.
Do not merge into `main` or claim detection significance. No upper limit
calculation has been added.

## One-command native C PS -> vendor likelihood -> candidate smoke workflow

Added `tools/run_ape_synthetic_candidates.py`, `src/ape_c_lpsd/spectral/short_synthetic_ps.py`, `src/ape_c_lpsd/analysis/synthetic_smoke.py`, `tools/stage_ape_ps_binary.sh` and `tests/test_ape_synthetic_run.py`. These provide a **synthetic-only, local Linux execution route**, distinct from the frozen historical 15,435-second / 5,035,956-bin APE profile.

### Stage original validated archive, without modifying it

The original frozen VIPER tarball is expected on the user's Mac at
`~/Downloads/APE_PS_VIPER_VALIDATED_20261022.tar.gz`, with archive SHA-256:
`0297b50652b06be90257e7471dc2cfce84cd15ed45add3129fdd34e0c2d6d48e`.

```bash
bash tools/stage_ape_ps_binary.sh \
  "$HOME/Downloads/APE_PS_VIPER_VALIDATED_20261022.tar.gz" \
  "$HOME/APE_PS_BINARY_STAGE_20261023"
```

The script prints the **different** executable SHA-256. The archived
`build/lpsd-exec` was compiled on Linux; a Mac may inspect and stage it,
but cannot be assumed to run it. Stage and verify on compatible Linux
(VIPER). Supply the executable's own SHA-256 to the JSON configuration.

### Run on a compatible Linux node with original vendor installed

Copy `configs/ape_synthetic_candidate_smoke.example.json` to a run-specific
JSON configuration. Provide a **new** `ape-synthetic-ps-*` run path, verified
`ps_executable`, its `ps_sha256`, unchanged `external_vendor_root`, a
physically justified QWP parameter configuration, and candidate index window.
The sample physics values are illustrative, **not user-approved APE constants**.

```bash
PYTHONPATH=src:. python tools/run_ape_synthetic_candidates.py \
  --config /absolute/path/to/synthetic_smoke.json
```

The requested explicit `synthetic_only=true`, `execute_c=true` controls
are mandatory. This path creates **new** calibrated-radian 60 s X/Y
coherent-tone noise fixtures; it does not use real data. For each quadrature,
the original native C executable computes one Method-3 FFT and three
Method-2 Constant-Q index partitions totaling **1301 bins**, with the same
numerical flags as the fixed original synthetic validation case. Seven-column
PS C outputs are split losslessly into legacy six-column PSD and two-column
PS; PS sidecars are validated and concatenated in index order. It then uses
the original vendor `fit_background.process_iteration` through
`fit_window` and original `stats.log_likelihood` through the APE adapter.
It scans a **caller-selected** set of candidate windows at a specified
finite coupling grid and writes `manifest.json` and `candidates.csv`.
`peak_shape_file=null` invokes vendor's original one-bin default `[1]`.
All PS used for the vendor model is rad², not rad²/Hz. Results are
**software-smoke candidates, not confirmed axion detections**.

**Known limitations:** This synthetic fixture is a coherent sinusoid, not
the physically broadened axion halo. Its arbitrary radian amplitude does
not correspond to an injected physical coupling, so it **cannot verify
coupling recovery**. User-supplied quadrature fractions and peak-to-RMS
normalization are not independently established. The vendor likelihood,
fit-quality and correlation/coverage assumptions have not been validated
with the actual C output. No external C binary or original vendor package
was available in the assistant's runtime; **no end-to-end native execution
can be claimed yet**. The committed focused tests also require execution
in the checked-out repository. Real-data processing, Slurm submission,
upper limits and physical false-alarm significance are deliberately excluded.

Focused checks (after checkout):
```bash
PYTHONPATH=src:. python -m unittest discover -s tests -p 'test_ape_synthetic_run.py'
PYTHONPATH=src:. python -m unittest discover -s tests -p 'test_ape_candidate_pipeline.py'
```
