# APE post-C-LPSD source audit (2026-10-31)

Scope: **only after original APE C-LPSD**. No upstream changes authorized.

Source of truth: user-provided original `Scalar-Dark-Matter-LPSD-1.0.0` files
`LPSDIO.py`, `lvk/fit_background.py`, `lvk/max_lkl_executable.py`,
`lvk/get_candidates.py`, `lvk/calc_upper_limit.py`, `lvk/stats.py`.

## Stage-by-stage trace

| Stage | Original vendor | Experimental APE | Difference |
|---|---|---|---|
| Native LPSD read | `LPSDOutput` | same original function, twice per observation | Two quadratures required |
| Combined observation | none | `logaddexp(logPS_X, logPS_Y)` | **APE required** |
| HDF5 serialization | `LPSDData.save_to_hdf5` | native HDF5 schema written with `h5py`, reread by original `LPSDOutput` | **Necessary vendor interface workaround** (below) |
| Background | `fit_background_in_file`, `LPSDJSONIO`, original skew-normal fitting | same functions, original CLI defaults | no statistical difference intended |
| Group loading | `LPSDDataGroup` | same | single detector label `APE` |
| Likelihood input | `make_args` | same | `stats.get_beta` APE response override |
| Likelihood/MLE | `process_q0` -> `calc_max_lkl` -> `stats.log_likelihood` | same | none intended |
| MLE serialization | `MaxLKLIO` | same | none |
| Candidates | `get_candidates.main` | original AST prefix up to `global_mngr`, including `cluster` | omit H1/L1-specific suffix |
| UL numeric | `get_upper_limits_approx` | same per frequency | vendor `main` lookup-table acceleration, smoothing, LIGO-specific comparisons and strain conversion **not reused** |

## Verified vendor interface issue

`LPSDIO.LPSDData.__init__(logPSD, freq=..., **kwargs)` assigns
`self.freq` only if `freq is None`. Supplying the original native frequency
array and then calling `save_to_hdf5` raises
`AttributeError: 'LPSDData' object has no attribute 'freq'`.
This was reproduced against the uploaded original source in a local isolated
vendor package. The combined-output adapter therefore writes the **same vendor
HDF5 schema** (`frequency`, `logPSD`, original attributes), preserving NaN
frequency gaps and metadata; then rereads using original `LPSDOutput`.

Original `fit_background.parse_args` supplies required `**kwargs`
(`kernel_size`, `verbose`, `plot_path`, `bin_scaling`,
`cfd_fraction`, `k_min`, `k_pruning`, `k_plateau`).
Calling `fit_background_in_file` with no options fails in
`process_iteration`; unified runner now supplies original CLI defaults.

## Critical scientific and completion checks NOT passed

1. `ape_beta_grid_npz` provenance/units and power-response normalization have
   not been numerically verified against the validated APE optical model.
2. Original C-LPSD output pairs and an APE production sample are unavailable
   in this runtime; no end-to-end post-LPSD fit/MLE/candidate/UL run completed.
3. `findiff` is missing from the local runtime, preventing direct import of
   original `lvk.stats` and candidate integration testing.
4. Original vendor UL `main` builds an interpolated 1000-point lookup table,
   then applies smoothing and LIGO-specific plotting/conversion. The unified
   APE driver calls its unchanged `get_upper_limits_approx` directly instead.
   This is a remaining execution-method difference, even though the underlying
   numerical function is original. Do **not** claim exact workflow parity.
5. Candidate extraction uses an AST-truncated original `get_candidates.main`.
   This preserves original statistical statements but is not a direct original
   entry-point call. It is brittle and requires a tested single-detector entry.
6. The original `cluster` implementation has a trailing-cluster behavior
   inherited unchanged. It must **not** be silently patched in a fidelity
   project.
7. This branch still contains superseded experimental adapters. They were
   not deleted without proving that no other experimental entry points depend
   on them.
8. A sum of two channel powers changes the null distribution; fitting the
   original vendor skew-normal model to that sum is the requested method,
   but empirical calibration/coverage has not been established.

## Status

**NOT CERTIFIED** as the closest possible implementation. The source-level
audit found real mismatches and corrected two interfaces, but end-to-end
integration and strict minimal-difference certification remain outstanding.
Do not use experimental outputs as scientific upper limits.
