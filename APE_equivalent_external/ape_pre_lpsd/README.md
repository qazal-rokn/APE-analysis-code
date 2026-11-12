# APE pre-LPSD stages (verbatim import of existing APE implementation)

This directory holds the **APE-specific inputs and preparation stages** that
precede the original external C-LPSD estimator. All copied Python modules,
command-line tools, tests, configs and historical documentation retain the
**exact Git blobs (unchanged bytes and executable modes)** from
`qazal-rokn/APE-analysis-code` at commit
`03722edb975e91a3653d3da673e5a521bb3a2854` (the source `main` snapshot).

This is **not** a rewrite of the external program, and it introduces no new
mathematical or statistical algorithm. The original vendor C/Python source
remains unchanged in the sibling `APE_equivalent_external/LPSD/` and
`APE_equivalent_external/scalar_dark_matter_lpsd/`, and the immutable reference
still lives at `vendor_direct_copy/original_program/`.

## Scope / boundaries

1. Read exactly five synchronized Dewesoft channels (X/Y at 1f and 2f,
   and cavity transmission) without changing sample timing or values.
2. Apply the **already established** APE quality cuts and physically
   continuous-segment selection. Preserve original indices, never fill gaps.
3. Apply the **existing** APE RMS calibration to produce separate
   `rho_x` and `rho_y` angle series, plus `i0_w` and `eta0`.
4. Write the frozen, calibrated segment as float64 HDF5 with the original
   provenance and non-overwrite rules. Do **not** combine X/Y in time.
5. Stop at the **input boundary to the original C-LPSD**. No C changes, spectral
   renormalization, background model, likelihood, or upper limits here.

Key entrypoints copied without editing:

- `tools/ingest_dxd.py`: validate a raw recording.
- `tools/inspect_cleaning*.py` and `tools/rank_cavity*.py`: diagnostics.
- `tools/find_longest_good_segment.py`: read-only candidate search/report
  (does **not** itself generate the frozen manifest).
- `tools/run_production_automatic_review.py`: older, review-only path that
  can produce a frozen manifest; requires the configured external reference.
- `tools/validate_frozen_segment_calibration.py`: inspect calibration of a
  short slice of an approved frozen manifest.
- `tools/calibrate_frozen_segment_to_hdf5.py`: create the calibrated HDF5
  product from an approved frozen manifest plus its SHA-256 provenance.

## Output and handoff

```text
APE .dxd
  -> read-only five-channel ingestion
  -> cavity/PEM cleaning + verified continuous segment
  -> RMS calibration (X/Y independent)
  -> canonical float64 HDF5: i0_w, eta0, rho_x, rho_y
  -> original C-LPSD executable (separate later execution boundary)
```

No execution or integration into the C-LPSD entrypoint has been performed as
part of this import. Existing real data and historical outputs are untouched.

## Tests and environment

From this directory, in a configured Python environment:

```bash
python3 -m pip install numpy PyYAML h5py matplotlib
PYTHONPATH=src python3 -m unittest discover -s tests -p 'test_*.py'
```

These are **imported existing synthetic tests**, not proof they passed in the
new directory. They require execution in a suitable environment. Real DXD
I/O additionally requires the Dewesoft Python header and native library,
located through the configured `reference_project`. Some legacy/reference
review tools also import `ape_analysis` from a separate original APE project;
that external project has **not** been copied here. Do not mistake an available
adapter for a self-contained reference implementation.

The configuration files are source examples and historical settings, not
universal up-to-date instrument calibrations. In particular, verify the
run-specific `cavity_transmission_w_per_v`, photodiode response and channel
names before any real-data work. The historical automatic-review YAML contains
machine-specific paths; use a separate run-local config rather than editing
this archived copy.

## Fidelity checks

See `IMPORT_MANIFEST.json` for every tracked file and its original Git blob
SHA. Identical Git blobs prove source-file identity, **not** successful
execution, calibration correctness for new data, or end-to-end equivalence.

The external post-LPSD stages will be analyzed and changed separately.
