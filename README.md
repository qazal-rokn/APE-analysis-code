# APE C-LPSD Comparison Project

This is a separate, approval-gated project for comparing the established APE
Python/SpiCyPy LPSD workflow with the compiled C `lpsd-exec` backend.

Current status: **contract proposed; no scientific implementation exists**.
Controlled synthetic cache/transport validation is implemented; real-data runs
remain refused. See `docs/QUICK_START.md` and `docs/REPRODUCIBILITY.md`.

Run the compiled-C synthetic validation with
` .venv/bin/python tools/run_synthetic_validation.py `. It uses temporary
synthetic HDF5 only and prints an auditable report.

For persistent dense-grid synthetic artifacts and measured small-grid resource
benchmarks, run `.venv/bin/python tools/run_dense_synthetic_benchmark.py`.

The intended invariant is that existing APE loading, channel mapping,
calibration, quality cleaning, and valid-data definitions remain unchanged. The
only intended replacement is the time-domain-to-frequency-domain estimator.

Authoritative proposal: `docs/ANALYSIS_CONTRACT_C_LPSD_BACKEND.md`.
Implementation must wait for explicit researcher approval and resolution of the
open decisions in that contract. A later, separate approval is required before
processing real experimental data.

Reference inputs are read-only:

- `/Users/user/Desktop/Code/Projects/APE-data-analysis`
- `/Users/user/APE-main-code-inbox`
- `/Users/user/Downloads/Scalar-Dark-Matter-LPSD-1.0.0/LPSD`
- `/Users/user/APE_Local_Data/After_wedged_EOM/AXION_COHERENCE_CLEANED_LPSD/CALIBRATED_CACHE`

The project must never write into those locations.

## Proposed layout

```text
.codex/agents/     governed specialist roles
configs/           approved run and backend settings
docs/              contract, method documentation, run records, changelog
notebooks/         optional presentation notebooks; no authoritative logic
results/           generated comparison artifacts; source data never stored here
src/               small import, segment, bridge, parser, plot, provenance modules
tests/             synthetic and format-validation tests
tools/             thin, auditable command-line entry points
```

For the approved C-native 65–9000 Hz planning stage, use `tools/plan_c_native_grid.py`. It reads only configuration and writes an immutable no-data grid manifest; it never exports signal data or invokes C. See `docs/C_NATIVE_GRID_PLANNER.md`.

## Planned low-frequency FFT extension

A separate future module, `src/ape_c_lpsd/fft_extension/`, will provide a
one-channel-at-a-time, full-uninterrupted-record, one-sided voltage FFT PSD for
low-frequency coverage. It does not replace the C-LPSD backend and does not
claim numerical identity with it.

```text
validated X voltage ──> FFT extension ──> separate fft_full_record PSD/ASD
                    └─> C-LPSD backend ─> separate c_lpsd PSD/ASD
validated Y voltage ──> FFT extension ──> separate fft_full_record PSD/ASD
                    └─> C-LPSD backend ─> separate c_lpsd PSD/ASD

separate products ──> optional display-only combined plot
```

The future FFT's display range, native FFT-calculable range, and C-LPSD
calculable range are separate configured concepts. See
`docs/FFT_EXTENSION.md`, `docs/FFT_EXTENSION_ARCHITECTURE.md`, and
`docs/FFT_EXTENSION_VALIDATION_PLAN.md`. The implementation and its tests have
not been run on real data or as a validation suite.

## Ingestion and longest-segment selection

The independent `src/ape_c_lpsd/ingestion/` package supports a future manual
mode (explicit ordered paths/slices) and automatic mode (an explicit
chronology-evidenced source manifest). It reuses the frozen reference quality
and calibration functions, creates review artifacts, splits at every invalid
sample or unverified boundary, and freezes the globally longest valid segment
for later separate X/Y FFT or C-LPSD products. It never infers physical order
from filenames or modification times. See
`docs/INGESTION_AND_SEGMENT_SELECTION.md` and `docs/INGESTION_ARCHITECTURE.md`.

## Calibrated quadrature review

The review-only `calibrated_review` stage applies the frozen reference scalar
calibration denominator independently to the two recorded voltage quadratures:
`X_ellipticity = pem1_x / D` and `Y_ellipticity = pem1_Y / D`.  It preserves
the exact `raw_index // 1000` reduced-record mapping and retains invalid
locations as gaps.  It makes no phase, X/Y gain-imbalance, or leakage
correction, and it does not combine X/Y.  Its outputs are dimensionless
ellipticity; `D` is recorded in V.  These review plots are not spectra and do
not authorize HDF5 export or C-LPSD execution.
