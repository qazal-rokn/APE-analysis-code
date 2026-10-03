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
