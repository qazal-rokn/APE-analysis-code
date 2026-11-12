# Automatic Ingestion, Cleaning Review, and Segment Selection

## Modes

- **Manual:** researcher declares ordered paths and exact raw slices; the package
  freezes them as a `frozen_validated_segment` manifest.
- **Automatic:** requires a researcher-provided/explicit source manifest with
  chronology evidence for every source. Filename order and filesystem time are
  rejected as physical-continuity evidence.

## Reused reference workflow

The adapter loads from the configured read-only reference project and calls,
without replacement:

1. `ape_analysis.quality.cavity_lock.build_cavity_lock_mask` with
   `CavityLockParameters`;
2. `ape_analysis.quality.pem_modulation.build_pem_modulation_quality` with
   `PemModulationParameters`;
3. `ape_analysis.calibration.polarization.build_polarization_calibration` with
   `PolarizationCalibrationParameters`.

Default YAML settings are the exact current reference defaults: cavity lock
`f_lock=.90`, `60 s`, quantile `.90`, three refinements; PEM settings
`20 mW/V`, `.7`, `.7 A/W`, `1e6 V/A`, `sqrt(2)`, `eta0>=.019`; and
polarization settings `20 kHz`, `.05 s`, `1000` raw/reduced mapping, `.7 A/W`,
`1e6 V/A`, `sqrt(2)`. The adapter refuses settings that differ from those
frozen defaults. A future approved threshold change must use a new YAML file
and is recorded in its settings snapshot/provenance; it cannot overwrite an
existing review run.

## Continuity rule

The automatic selector splits at every invalid raw sample. It joins two runs
across files only when all are true: the first reaches its file endpoint, the
next begins at raw zero, the files are adjacent in the explicit manifest, and
the incoming boundary is explicitly verified. Missing files, sample-rate
mismatch, overlaps, source gaps, or unverified boundaries terminate a segment.

## Planned review products

Each non-overwriting review run stores source manifest/checksums, configuration
snapshot, metadata and per-file fractions, mask/segment summaries, selected
segment YAML/JSON, logs, and a human-readable report. Future raw and cleaned
plots use labelled min/max display envelopes; decimation is visual-only and is
never supplied to cleaning or spectra. Diagnostic channels are user-configured,
including cavity, PEM harmonics/quality, and X/Y.

The C backend gate accepts only a `frozen_validated_segment` manifest from one
of these modes. It does not accept a loose path or a guessed segment.

## Implemented plotting products

`ingestion.plotting.generate_ingestion_review_plots` creates one UUID-named
review directory. For every configured diagnostic channel it writes separate
`raw/`, `cleaned/`, and `selected_longest_segment/` directories with a
non-overwriting overview PNG and metadata. Cleaned stage adds a valid-fraction
mask plot. Selected stage may add a bounded full-resolution focus/zoom plot.
Overview rendering uses min/max envelopes labelled with their original-sample
width; this is display-only and never changes the supplied data, mask, or
segment. `plot_metadata.json` records output paths/hashes, source identifiers,
units, rates, stage, display method, and config snapshot.
