# APE vendor-faithful pipeline (experimental, NOT VALIDATED)

This branch is a full copy of `development/likelihood`, with experimental additions isolated from the existing validated pipeline. No original vendor files or production APE files are changed.

## Strict scientific contract

From the C-LPSD output onward, use original vendor functions and workflow wherever possible. The ONLY scientific replacement is the detector response (vendor LIGO `stats.get_beta(freq_Hz, f_A_star[ifo])` versus calibrated APE optical response). Data-format adapters may be necessary but must not change statistical meaning. Preserve original vendor behavior even when suspicious; log any issues for validation rather than silently correcting them.

## What has actually been added

`src/ape_c_lpsd/analysis/ape_vendor_faithful_experimental.py` delegates to the unmodified vendor `fit_background_in_file` for complete vendor segmentation, fitting and JSON serialization, provided a vendor-readable LPSDOutput input; and delegates to unmodified vendor `calc_max_lkl` for any normalized single- or multi-bin peak shape. These are building blocks, not a complete runnable end-to-end analysis. UNTESTED.

## Remaining work before the requested complete parallel pipeline exists

1. Convert native APE C-LPSD PS and available observation metadata into the vendor's `LPSDOutput` / `LPSDDataGroup` format without changing power units or frequency grid. Validate conversion parity and segment identity.
2. Reuse vendor `fit_background_in_file` with vendor defaults/settings, original segment boundaries and JSON metadata, rather than APE's existing 2,000-bin segment scheme. Verify fit outputs.
3. Reproduce vendor `max_lkl_executable.make_args` behavior exactly for multi-bin windows: skip NaN records, look up background by record metadata and tested frequency, apply `chi_sqr >= max_chi_sqr` cut, evaluate spline, assemble skew-normal args, preserve frequency bounds and peak-shape indexing. Replace ONLY `stats.get_beta(freq_Hz, f_A_star[ifo])` with APE's calibrated signal-power normalization. Ensure the physical meaning of mu is consistent.
4. Preserve original `calc_max_lkl` exactly, including signed-mu validity search, seed grid, Nelder-Mead and Fisher sigma. Note: vendor `get_valid_mu_ranges` omits explicit peak_shape; do not patch it here.
5. Preserve vendor result format and original `get_candidates.py` methods where scientifically applicable. Its H1/L1 consistency test is detector-specific and cannot be relabeled X/Y without a separate physical argument. This is a necessary detector-geometry exception, not evidence of identical observatories.
6. Preserve vendor `get_upper_limits_approx` and input conventions; replace physical coupling conversion only.
7. Run identical-input numerical parity, segment/quality-cut parity, synthetic injection, and end-to-end tests. No tests have been executed for this branch.

## Known non-equivalences in existing validated pipeline

- `tools/ape_real_full_pipeline.py` performs APE-defined 2,000-bin background segmentation and aborts on failed fits; vendor `fit_background_in_file` uses its own segmentation/metadata and vendor `make_args` skips failed records.
- Existing APE likelihood is single-bin only; vendor supports multi-bin.
- Existing APE candidate clustering and trial handling are not the original vendor's complete procedure.
- Existing APE treats two quadratures as records; vendor records represent time segments and detectors. Statistical independence and correlation assumptions are not established.

## Acceptance

Status: EXPERIMENTAL, INCOMPLETE, UNTESTED. This branch must not be presented as a completed vendor-identical pipeline, nor used for published physical exclusions. The existing APE pipeline is unchanged.
