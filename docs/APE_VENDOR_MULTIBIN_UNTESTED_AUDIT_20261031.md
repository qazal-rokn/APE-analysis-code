# Vendor multi-bin likelihood bridge: UNTESTED audit (2026-10-31)

Branch: `feature/vendor-multibin-likelihood-untested`. This branch intentionally leaves all existing production analysis, raw data, vendor code, and verified single-bin adapters unchanged.

## Added, not validated

`src/ape_c_lpsd/analysis/ape_vendor_multibin_untested.py` passes vendor-format `Y`, `bkg`, `model_args`, `peak_norm`, and `peak_shape` directly into the original `scalar_dark_matter_lpsd.lvk.max_lkl_executable.calc_max_lkl`. Vendor source hash is checked against the existing APE pinned SHA256. No vendor likelihood, optimizer, background, or statistics code is copied or altered.

This is a callable **isolated adapter**, NOT a connected end-to-end multi-bin pipeline. It has not been executed or tested.

## Source comparison

- Original vendor `max_lkl_executable.make_args`: extracts `logPSD[:, freq_idx:freq_idx + len(peak_shape)]`, evaluates the background spline at those frequencies, uses one skew-normal parameter triplet per record, computes `peak_norm` using vendor `stats.get_beta`, and passes all arrays to `calc_max_lkl`.
- Original vendor `stats.log_likelihood`: combines `exp(bkg) + peak_norm * peak_shape * mu_DM` and sums log skew-normal PDF contributions over records and bins.
- APE `ape_vendor_fit.py`: explicitly enforces `Y.shape[1] == 1` and `peak_shape == [1.0]`; it remains untouched.
- APE `vendor_statistics.fit_window`: calls original vendor `fit_background.process_iteration` and `bkg_model` unchanged, but fits an APE-selected local window, with APE-configured settings. It is not yet demonstrated to reproduce the vendor's full segmentation, quality filtering, or all runtime settings.
- APE `vendor_statistics.log_likelihood`: delegates to original vendor `stats.log_likelihood`; this does not by itself establish full-workflow parity.
- APE `vendor_statistics.factor_power_template`: constructs a shape/norm decomposition from an APE power template. Whether it matches the external program's actual `peak_shape` and `get_beta` convention for the intended APE response remains to be demonstrated.

## Known unverified or different integration behavior

1. The original vendor `make_args` constructs and filters multi-bin windows from its own data and background metadata; the new adapter expects prepared arrays and does NOT reproduce this preparation.
2. The original vendor full workflow includes its own background segmentation, skew-normal fitting, quality cuts, frequency iteration, peak-shape loading, and output. These are NOT wired into the new adapter.
3. The original vendor's `peak_norm` uses `get_beta(freq_Hz, f_A_star[ifo])`. APE apparatus response must replace the physical normalization only, with its precise units and `mu` definition documented.
4. Adjacent APE Constant-Q bins may be correlated. Numerical parity of the original likelihood is not statistical coverage validation.
5. No multi-bin parity, injected-signal recovery, local significance calibration, or observed-limit coverage tests have run.

## Required acceptance before use

Review vendor `fit_background.py`, `max_lkl_executable.py`, `stats.py`, peak-shape construction, `make_args`, and APE data preparation line by line. Run identical-input parity tests with vendor `calc_max_lkl`; verify background parameters and all selection/filter conventions; validate apparatus-only differences; then test end-to-end on controlled fixtures. Do not claim full scientific equivalence or production readiness before these checks.
