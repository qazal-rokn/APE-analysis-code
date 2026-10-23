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
