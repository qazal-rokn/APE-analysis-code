# Analysis Contract: Whole-Acquisition Longest Clean Segment to Compiled C LPSD

Date proposed: 2026-09-24

No scientific backend is implemented by this contract. No synthetic or real APE
spectrum is authorized by its current status.

## REQUESTED

Create a separate project at `/Users/user/Desktop/Code/Projects/APE-data-analysis-c-lpsd` that preserves
the existing APE workflow through loading, channel mapping, calibration, quality
cleaning, and continuous-good-segment selection, but replaces only the Python
SpiCyPy/LPSD spectral estimator with the already-built compiled C
`lpsd-exec`. Search the complete chronological cached acquisition for every
maximal uninterrupted valid run, select exactly the globally longest run, export
only that run to temporary HDF5, invoke the C executable, parse frequency and
PSD, and create provenance-rich spectrum artifacts. Synthetic validation must
pass before a separately approved real-data run.

## INSPECTED CURRENT BEHAVIOR

### Reference APE project

The reference project currently:

- resolves Dewesoft channels and calibration through separate modules;
- forms approved quality masks and maps them to full-rate calibrated samples;
- defines a contiguous good segment as adjacent integer raw indices within one
  source identity;
- never fills or bridges index gaps;
- has an approved Python native-grid, one-prefix LPSD adaptation using a Kaiser
  window and one-sided normalization `2/(fs * sum(w^2))`; and
- documents a prior first-file segment of 68,817,000 samples at 20,000 Hz, but
  that result is not a whole-cache global-longest determination.

This new project will not import that prior spectrum or assume its segment is
globally longest.

### Inspected calibrated cache

`cache_metadata.json` identifies cache version 2, float64 samples, 20,000 Hz,
246 nominal 600-second blocks, 2,952,000,000 positions over the selected time
range, 215,999,999 observed samples, and 1,220 bad intervals. Each block names
`rho_x`, `rho_y`, and `valid` `.npy` files and records global indices, absolute
times, observed sample count, and source-derived subsegments. The source list is
chronological and omits file suffix `_0002`; later blocks may have no observed
samples. Therefore consecutive block IDs are not sufficient evidence of physical
continuity. Continuity must be proven from metadata/index/time relationships and
the valid mask at every boundary.

The cache's calibration signature records constants and a cavity-mask hash. The
new project will consume these as provenance, not recompute or alter them.

### Inspected C LPSD program

The inspected executable is:

```text
/Users/user/Downloads/Scalar-Dark-Matter-LPSD-1.0.0/LPSD/manual-build/lpsd-exec
SHA-256 fc8a912fba75aee0f26806658ace5414b01ac8887c6fe2675f657b0bf10030ef
Mach-O 64-bit x86_64
embedded LPSD_VERSION: LPSD 28.06.2007
```

The inspected source tree describes a 2023 scalability adaptation and CMake
project version 0.0.1. It links HDF5, GSL, and CBLAS/OpenBLAS. CMake builds
`lpsd-exec` from the listed C sources and links `HDF5::HDF5`, `GSL::gsl`, and
`GSL::gslcblas`. A reproducibility record will pin the executable and relevant
source-file hashes and capture `otool -L`/platform details.

The executable supports a named HDF5 dataset with `-D` (source default
`strain`), sampling rate `-f`, LPSD mode `-h 0`, input/output `-i/-o`, frequency
limits `-s/-t`, total/requested bin controls `-J/-n/-N`, window `-w`, peak
sidelobe level `-p`, overlap `-l`, detrending/linear regression `-r`, desired and
minimum averages `-a/-m`, input scale `-x`, and memory control `-M`.

In LPSD mode, the source opens an HDF5 dataset as native doubles and writes a
text result. Output columns are selected by the configured `FORMAT`; `f` is
frequency, `D` PSD, `d` square-root PSD, `N` averages, `R` resolution bandwidth,
and `b` bin number. The implementation computes one-sided PSD normalization in
the C estimator; the Python bridge must not add a second normalization. Effective
configuration cannot be left to an ambient `LPSDCFN` or interactive defaults.

## UNDERSTANDING

The scientific invariant is:

```text
approved cached calibrated samples + approved cached valid mask
  -> complete chronological continuity audit
  -> exactly one globally longest maximal valid, physically continuous run
  -> lossless float64 HDF5 transport
  -> compiled C lpsd-exec only
  -> parsed native frequency and PSD
  -> ASD only as sqrt(PSD) for display/output
```

No calibration occurs in the C program. Input quantity and unit propagate
unchanged into the spectrum: if `x` has unit `U`, PSD has `U^2/Hz` and ASD has
`U/sqrt(Hz)`.

## PLANNED IMPLEMENTATION AFTER APPROVAL

### 1. Read-only cache adapter

Validate metadata schema, cache signature, expected file set, array shape/dtype,
block/global-index consistency, sample rate, and read-only paths. Stream arrays
with memory mapping; never rewrite cache files. Preserve chronological order.

Inputs: `cache_metadata.json` and named block arrays.

Outputs: validated block descriptors and streamed float64 signal/mask views.

### 2. Whole-acquisition continuity inventory

Scan every block in chronological global-index order. A run continues across a
block or source-file boundary only if all of the following are true:

1. both boundary samples exist and are valid;
2. global sample indices differ by exactly one;
3. reconstructed timestamps differ by exactly `1/fs` within an approved
   floating-point tolerance derived from metadata precision;
4. source-derived metadata contains no missing interval at the boundary; and
5. the sample rate and signal definition are unchanged.

Any failed condition closes the run. Missing observations, absent source files,
invalid-mask values, index jumps, time jumps, metadata inconsistencies, or
nonfinite selected signal values are boundaries, never values to repair. Emit a
complete segment inventory with endpoints, sample count, duration, contributing
blocks/source files, boundary evidence, and selection rank. Select the unique
largest sample count; use the approved tie rule if counts tie.

### 3. Selected-segment stream

Read exactly the selected half-open global interval in original order and assert
that every mask value is true and all continuity checks still pass. Do not join,
fill, interpolate, resample, smooth, detrend, or otherwise transform it in the
Python layer. Any C-side detrending/windowing is part of the explicitly approved
estimator configuration and will be recorded.

### 4. Temporary HDF5 export

Create a run-scoped temporary HDF5 file, not a source/cache artifact. Write one
rank-1, float64 dataset under an approved explicit dataset name, chunked only as
an I/O representation with no numerical transformation. Store transport
attributes (quantity name, unit, sample rate, count, segment endpoints, hashes)
for the wrapper's checks; the C program reads only the numeric dataset. Verify
shape, dtype, order, endpoints, and a streaming content hash after writing.

### 5. Compiled C invocation

Invoke the pinned executable directly with an argument vector (no shell), a
controlled environment, an explicit project-owned LPSD configuration, explicit
input dataset, sample rate, method, output path, window/PSLL, overlap,
detrending, frequency limits/grid/bin controls, averaging, scaling, memory, and
batch iteration. Capture exact arguments, executable/source/config hashes,
stdout, stderr, exit status, runtime, platform, and dynamic-library provenance.
Reject nonzero exits or unexpected files. The command must genuinely execute
`lpsd-exec`; no Python estimator fallback is permitted.

### 6. Output parsing and normalization guard

Use a project-owned configuration whose `FORMAT` starts with explicit frequency
and PSD columns and includes diagnostic averages/resolution/bin fields if
approved. Parse headers and finite numeric rows strictly; verify row count,
monotonic frequency, requested range, nonnegative PSD, and diagnostic fields.
Preserve C PSD values without renormalization. Compute ASD solely as
`sqrt(PSD)` and label units from the input quantity.

### 7. Artifacts, plot, and provenance

Write versioned, non-overwriting results containing native C frequency, PSD,
ASD, diagnostic columns, units, segment inventory/selection record, effective
configuration, exact command, hashes, logs, and run record. Plot PSD and/or ASD
without smoothing or interpolation. Notebooks may display outputs but will not
contain authoritative processing logic.

### 8. Modular boundaries

Planned modules are cache inspection, chronology/continuity, segment selection,
HDF5 transport, C process invocation, C output parsing, spectrum artifact I/O,
plotting, and provenance. `src/physics/` will contain only approved unit and
normalization equations needed to keep code and documentation synchronized.

## SCIENTIFIC DEFINITIONS AND UNITS

For selected samples `x_n` with unit `U` and sample rate `f_s`:

- no Python-layer calibration or scaling is applied;
- C input scale must be `-x 1` unless the researcher explicitly approves
  otherwise;
- C output `D` is interpreted as one-sided PSD `S_x(f)` in `U^2/Hz`, subject to
  synthetic normalization validation;
- displayed ASD is `sqrt(S_x(f))` in `U/sqrt(Hz)`;
- frequency is in Hz; and
- duration is `N/f_s`, while endpoint span is `(N-1)/f_s`.

The final equations for the approved C mode, window, detrending, averaging, and
grid will be transcribed from the inspected source into
`docs/EQUATION_REFERENCE.tex` and mapped to implementation after approval.

## SCIENTIFIC RISKS AND CONTROLS

- **False cross-gap stitching:** require index, time, observation, mask, and
  metadata boundary evidence; test missing-file and block-boundary cases.
- **Wrong cached signal:** require researcher approval of `rho_x`, `rho_y`, or
  another exact definition before implementation.
- **Cache not identical to reference modules:** pin cache signature/constants and
  document that the project consumes, rather than recomputes, cached calibration
  and masks; verify provenance against the approved source.
- **Ambient C configuration:** use an explicit hashed config and controlled
  `LPSDCFN`; record effective settings printed by the program.
- **C output-column ambiguity:** set and test explicit `FORMAT`; parse by verified
  header, never positional guesswork alone.
- **Normalization mismatch/double scaling:** validate tone and noise analytically,
  keep `-x 1`, and never renormalize C PSD in Python.
- **Window/detrending/averaging changes:** require explicit values and record C
  effective settings; tests verify their effects.
- **Frequency-grid batching/checkpoint behavior:** test small multi-batch jobs,
  completeness, ordering, duplicate/missing bins, and clean run-scoped outputs.
- **Executable/source drift:** hash binary, source, config, dependencies, and
  invocation for every run.
- **Non-determinism or stale checkpoint reuse:** fresh unique work directory and
  output names; byte/numeric deterministic-rerun tests.
- **Architecture/dependency portability:** document x86_64 binary and linked
  library versions; provide reproducible CMake build instructions without
  rebuilding or installing dependencies unless separately approved.

## WILL NOT DO

- No modifications to either reference project, notebooks, C source/executable,
  source data, cache files, or existing APE results.
- No raw `.dxd` reprocessing in the proposed implementation.
- No recalibration or redefinition of the cached valid mask.
- No gap stitching, fill, interpolation, resampling, smoothing, filtering,
  normalization, or unapproved rejection.
- No Python/SpiCyPy LPSD result copied, relabeled, or used as the production
  backend.
- No fallback estimator.
- No candidate finding, vetoes, axion limits, upper limits, coupling conversion,
  or physics interpretation.
- No real-data estimator run until synthetic validation passes and the researcher
  gives a second explicit approval.
- No package installation or C rebuild under this contract.

## PROPOSED FILES FOR THE APPROVED IMPLEMENTATION

Scaffold/governance files already created by the first deliverable are not
scientific implementation. After approval, expected additions/changes are:

```text
configs/cache.yaml
configs/c_lpsd.yaml
src/ape_c_lpsd/cache.py
src/ape_c_lpsd/segments.py
src/ape_c_lpsd/hdf5_transport.py
src/ape_c_lpsd/c_backend.py
src/ape_c_lpsd/output.py
src/ape_c_lpsd/plotting.py
src/ape_c_lpsd/provenance.py
src/physics/spectral_equations.py
tests/test_cache_format.py
tests/test_segment_selection.py
tests/test_hdf5_transport.py
tests/test_c_backend_synthetic.py
tests/test_output_units.py
tests/test_determinism.py
tests/test_data_order.py
tools/run_synthetic_validation.py
tools/run_c_lpsd.py
docs/CACHE_INPUT_FORMAT.md
docs/WHOLE_ACQUISITION_SEGMENTS.md
docs/HDF5_EXPORT.md
docs/C_LPSD_INVOCATION.md
docs/C_LPSD_OUTPUT.md
docs/UNITS_AND_NORMALIZATION.md
docs/PLOTTING_AND_ARTIFACTS.md
docs/PROVENANCE.md
docs/VALIDATION.md
docs/LIMITATIONS.md
docs/REPRODUCIBILITY.md
docs/EQUATION_REFERENCE.tex
docs/SCIENTIFIC_CHANGELOG.md
README.md
```

No file in `/Users/user/Desktop/Code/Projects/APE-data-analysis` would change.

## CONTROLLED SYNTHETIC VALIDATION PLAN

Before any real-data estimator run:

1. **Tone recovery:** generate a deterministic on-grid and documented off-grid
   sinusoid with known amplitude/unit; verify peak frequency and C PSD/ASD
   normalization using the approved window and detrending settings.
2. **Noise-only behavior:** seeded white Gaussian noise with known one-sided PSD;
   test ensemble/band statistics with tolerances justified by the number of
   averages and window, not exact-bin equality.
3. **HDF5 round trip:** verify rank-1 float64 dataset name, shape, order, endpoint
   values, and streaming hash; confirm the C executable reads it.
4. **No-gap stitching:** masks and metadata containing one-sample, multi-sample,
   index, timestamp, missing-observation, and missing-file gaps must always yield
   separate runs.
5. **Cross-file global longest:** sequential synthetic files/blocks where the
   longest run is not in the first file, plus a valid run that legitimately
   crosses a proven continuous boundary.
6. **Tie behavior:** exercise and verify the researcher-approved equal-length
   tie rule.
7. **Output units:** voltage-labelled input yields `V^2/Hz` and `V/sqrt(Hz)`;
   calibrated polarization/ellipticity input yields its exact approved
   unit-squared/Hz and unit/sqrt(Hz), with no numerical unit conversion.
8. **Deterministic reruns:** identical input/config/executable produces identical
   frequency/PSD values and stable scientific hashes; timestamps/runtime may
   differ only in provenance fields.
9. **Data-order preservation:** monotonic sentinel sequences spanning blocks and
   boundaries arrive in HDF5 exactly once and in original order.
10. **Output parser:** reject malformed headers, wrong formats, NaN/negative PSD,
    nonmonotonic or duplicate frequencies, missing batches, and unexpected rows.
11. **Command/provenance:** verify direct invocation, complete explicit options,
    controlled environment, captured logs/status, and pinned hashes.
12. **Batch equivalence:** where batching is required, compare a small one-job
    and multi-job synthetic run for identical ordered native bins and PSD.

Passing criteria and numerical tolerances will be documented before tests run.
Synthetic test artifacts only may be generated after contract approval. Passing
tests authorize reporting validation results, not a real-data run.

## UNRESOLVED DECISIONS REQUIRING RESEARCHER APPROVAL

1. **Signal channel/definition:** export `rho_x`, `rho_y`, both as separate C
   runs, or another exact cached quantity? A magnitude must not be invented.
2. **Physical unit label:** confirm whether cached `rho_x/rho_y` are radians,
   dimensionless ellipticity, or another convention, and the exact display text.
3. **Tie rule:** if multiple maximal runs have equal sample count, select earliest
   chronologically (proposed) or stop for manual selection.
4. **Boundary tolerance:** approve a time-continuity tolerance. Proposed: compare
   integer global indices as primary truth and require reconstructed time error no
   larger than the cache metadata's documented reconstruction bound plus a small
   floating-point margin; never tolerate a missing sample.
5. **HDF5 dataset name:** proposed explicit `/strain` for compatibility with the
   C default, recorded with quantity/unit attributes.
6. **C window and PSLL:** exact `-w` and `-p` values. The C source defaults to
   Kaiser `-2`, PSLL 120 dB; the prior APE Python adaptation used 200 dB, while
   the C example uses 238.13 dB. These are not interchangeable.
7. **Detrending:** exact `-r` value: none (0), linear regression (1, C default),
   or endpoint line (2).
8. **Overlap and averaging:** exact `-l`, `-a`, and `-m`. The C defaults imply
   automatic overlap, 100 desired averages, and 10 minimum averages, which do not
   match the prior Python one-prefix/single-window policy.
9. **Frequency grid/range:** exact `fmin`, `fmax`, total bins `Jdes`, batch size,
   and whether to reproduce the prior APE coherence-based native-grid density or
   use the C program's approved native LPSD planning. The estimator comparison is
   not scientifically defined until this is fixed.
10. **Output format:** proposed `fDNRb` (frequency, PSD, averages, resolution,
    bin) in a project-owned config; approve diagnostic columns.
11. **C input scaling:** proposed `-x 1` with calibrated values already supplied.
12. **Batching/memory:** approve `-M`, `-n`, and `-N` strategy after a synthetic
    feasibility benchmark; these affect execution and possibly edge behavior.
13. **Comparison artifact scope:** whether this project should initially emit
    only the C spectrum or also ingest the existing Python result later for a
    side-by-side plot. The latter is not included in the first implementation.

## REAL-DATA GATE

A real cache scan may inspect metadata/schema only. Exporting real selected
samples to HDF5 or invoking `lpsd-exec` on them is prohibited until:

1. this contract and every method-defining unresolved decision are approved;
2. the implementation and documentation are complete;
3. every controlled synthetic validation passes;
4. physics and signal/statistics reviews report the result;
5. the final audit compares request, contract, code, validation, and docs; and
6. the researcher explicitly approves the real-data run in a separate decision.

## STATUS

**APPROVED FOR IMPLEMENTATION AND CONTROLLED SYNTHETIC VALIDATION ONLY.**

The approved profile uses cached `rho_x` and `rho_y` separately, the same
globally longest continuous valid segment, C Kaiser `-w -2`/`-p 200`, linear
regression `-r 1`, 50% overlap, `-a 100`, `-m 1`, `-x 1`, `fDNRb`, and metadata
derived $f_{min}=f_s/N$, $f_{max}<f_s/2$, and $J_{des}$ at $Q=10^6$. No real
experimental data may be exported or analyzed without a separate approval.
