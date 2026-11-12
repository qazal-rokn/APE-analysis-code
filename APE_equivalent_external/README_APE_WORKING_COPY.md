# APE equivalent of original Scalar-Dark-Matter-LPSD

This directory began as a **complete 68-file byte-for-byte original upstream
program** (upstream commit `f9e78cbe1031406e217c0bfb63a2f6f0043a0a5a`,
tree `237b3ba91539ca0ac188cc78c2906fb035b6c301`).

- `LPSD/`, `scalar_dark_matter_lpsd/`, `transfer_functions_O4/`,
  `examples/`, and all original root source/metadata files remain **unedited**.
- `ape_pre_lpsd/` is a **separate additional folder**, verbatim-importing the
  existing APE ingestion, cleaning, calibration, and calibrated-HDF5 interface.
- The immutable reference remains at
  `vendor_direct_copy/original_program/` outside this directory.
- Post-LPSD APE adaptations have **not** been applied.

Because this working copy contains an additional APE folder, its **overall**
Git tree SHA no longer equals the original root tree. This is expected.
Every original upstream file must still retain its original Git blob SHA;
use the immutable reference for exact numerical and source comparisons.

See `ape_pre_lpsd/README.md` and `ape_pre_lpsd/IMPORT_MANIFEST.json`.

**Do not** silently alter the vendor Constant-Q estimator, background fit,
skew-normal likelihood, maximization, candidate significance or upper limits.
All later APE-specific deviations require an explicit scientific justification,
a source diff and quantitative validation.
