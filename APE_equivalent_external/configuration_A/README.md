# Configuration A: ORIGINAL vendor single-bin likelihood for APE

**Status: IMPLEMENTED AS A SOURCE-PRESERVING ADAPTER; PHYSICAL RESPONSE NOT YET CERTIFIED.**

This configuration uses the original unmodified vendor Python code in
`../scalar_dark_matter_lpsd/` and original C LPSD in `../LPSD/`.
It uses the original `LPSDOutput`, `LPSDJSONIO`, `bkg_model`,
`calc_max_lkl`, `stats.log_likelihood`, `MaxLKLIO` and the original
single-bin `peak_shape = np.ones(1)`. The adapter replaces ONLY the
experiment-specific `peak_norm` (formerly LIGO `get_beta`) with an
externally supplied APE PSD response per unit fitted `mu`.

The pre-LPSD APE ingestion/calibration stage and PSD handoff remain in
`../ape_pre_lpsd/` and `../ape_post_lpsd/prepare_vendor_lpsd.py`;
we do not rerun or modify those accepted stages here. The background fit
is the unmodified vendor `fit_background.py` using existing accepted
X and Y background JSON files. This entry point fits **X and Y separately**:
joint X/Y statistical independence has not been established.

## Execute from APE_equivalent_external/

```bash
python -m configuration_A.run \
  --spectrum /absolute/path/APE_0_15435_X.h5 \
  --background-json /absolute/path/background_X_FULL.json \
  --response-npz /absolute/path/approved_ape_response.npz \
  --output-npy /absolute/path/config_A_X.npy \
  --channel X --ana-fmin 70 --ana-fmax 7000 \
  --min-log10mu -8 --max-log10mu 8
```

The `mu` scan exponents shown here are **illustrative only** and must be
selected from the verified units and normalization of the response NPZ.
For Y, supply a Y spectrum/background and `--channel Y`, with a distinct
output path. Never combine X/Y estimates without a validated covariance model.

The response NPZ contract is enforced by `ape_post_lpsd.ape_signal_response`:
`frequency_hz`, `peak_norm_x`, `peak_norm_y`, `metadata_json`.
Its exact frequency grid and PSD units must match C-LPSD. Its `mu_definition`,
response provenance, and optical response calibration require independent
physical validation. Labels in JSON are NOT proof of calibration.

## Hard boundaries

- **DO NOT** import any `ape_axion_halo_*`, `ape_native_kaiser_kernel`,
  `ape_provisional_qwp_multibin_vendor`, or `configuration_B` code.
- **DO NOT** put a multi-bin `peak_shape` into Configuration A.
- **DO NOT** change the vendor source; make any adaptation outside it.
- No official upper limit, sensitivity, or detection significance until APE
  injection/recovery, PSD normalization and statistical coverage checks pass.
- `calc_upper_limit.py` has LIGO-specific conversions/plot paths in its
  executable `main()`. The reusable `get_upper_limits_approx` function
  must be used unchanged with a scientifically validated `sigma_mu`; the
  physical coupling conversion belongs in a separate APE boundary.

Configuration B is documented separately in `../configuration_B/README.md`.
Its existing multi-bin/Monte-Carlo diagnostics are excluded from this
entry point and from Configuration A outputs.

For full source preservation, compare files under
`../scalar_dark_matter_lpsd` with the untouched vendor reference before
claiming byte-for-byte parity. This branch adds no changes inside that tree.
