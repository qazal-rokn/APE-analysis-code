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

## APE-only signal-power response generator (added)

The physical detector interface is now `configuration_A.build_response`.
It does **not** reimplement the vendor likelihood or the halo line shape.
It calculates the additive original-vendor `peak_norm` for each X/Y native bin
from the two **externally calibrated** instrument quantities:

- `R_k(f)`: measured PEAK rotation amplitude per unit axion-photon coupling,
  incorporating the two-QWP cavity response, PEM/lock-in calibration, and
  X/Y readout projection, in `rad / (GeV^-1)`;
- `K_k(f)`: actual native **C method-2 D output** in `Hz^-1` per
  unit time-domain angular variance, measured for the same bin and
  single-bin signal-power convention (not a generic uncalibrated Kaiser
  spectral-window function).

With `mu=(g/g_ref)^2`, the only new signal term is

`peak_norm_k(f) = 0.5 * [g_ref * R_k(f)]^2 * K_k(f)`

and `peak_shape=[1]` remains unchanged. The 0.5 applies because the
interface explicitly requests *peak* sinusoidal amplitude, not RMS amplitude.
Do not supply RMS or PSD-amplitude numbers as `R` without conversion.

The input `ape_calibration.npz` MUST contain the following five exact
array names, all on the full native C-LPSD grid:

```
frequency_hz
peak_rotation_rad_per_gev_inv_x
peak_rotation_rad_per_gev_inv_y
c_psd_hz_inv_per_unit_variance_x
c_psd_hz_inv_per_unit_variance_y
```

The accompanying `ape_calibration_provenance.json` must specify

```json
{
  "reference_coupling_gev_inv": 1e-10,
  "response_source": "REPLACE WITH ACTUAL MEASURED APE RESPONSE",
  "calibration_version": "REPLACE WITH ACTUAL VERSION",
  "native_c_kernel_source": "REPLACE WITH ACTUAL NATIVE C TEST",
  "single_bin_convention": "REPLACE WITH REVIEWED SOURCE CONVENTION",
  "peak_amplitude_calibration_source": "REPLACE WITH ACTUAL MEASUREMENT",
  "quadrature_projection_source": "REPLACE WITH ACTUAL CHANNEL CALIBRATION",
  "response_status": "PROVISIONAL_NOT_VALIDATED_FOR_PHYSICAL_LIMITS"
}
```

This JSON is a SCHEMA EXAMPLE, not a supplied physical measurement; do
not substitute placeholder text into a science run.

Run from `APE_equivalent_external/`:

```bash
python -m configuration_A.build_response \
  --calibration-npz /work/ape_calibration.npz \
  --provenance-json /work/ape_calibration_provenance.json \
  --output-npz /work/ape_A_response.npz
```

Pass the generated NPZ to `python -m configuration_A.run` with the
original C-LPSD spectrum and vendor background fit. The latter delegates
straight to the unmodified vendor `calc_max_lkl`. For post-fit
candidate screening and median-expected limit diagnostics, use existing
`ape_post_lpsd.ape_finalize_inference`, which calls original vendor
`cluster` and `get_upper_limits_approx`. It does not compute an observed
limit or perform LIGO's H1/L1 consistency veto.

The source contains **no completed numerical APE measurement or C-tone
calibration**. A correctly shaped NPZ is not a validated coupling response.
Until actual optical response and native PSD kernel measurements are
provided and injection/recovery and coverage are passed, outputs remain
CONDITIONAL diagnostics, not physical exclusion limits.

### Exact source-fidelity boundary

No changes have been made to `LPSD/`,
`scalar_dark_matter_lpsd/`, the original C numerical method, or
the vendor statistical functions. The APE-specific original data handoff,
background-fit handling, X/Y splitting, part combining, and inferred coupling
reporting were already implemented under `ape_post_lpsd/` and are reused.
The new builder changes *only* the external numerical `peak_norm` input.
Configuration B remains separate.


## Optional: construct the calibration from the EXISTING APE two-QWP model

`configuration_A.derive_qwp_calibration` imports the **unchanged existing**
`ape_post_lpsd.physics.qwp_signal.qwp_amplitude_per_coupling`, rather than
implementing another optical response equation. Its model remains
**PROVISIONAL** and is not equivalent to a measured complex cavity response.

Supply an `optics.json` containing exactly the original QWP model arguments
`finesse_s`, `tau_s`, `tau1_s`, `rho_dm_gev4`, `unit_conversion`,
with actual independently validated numerical values and units. Supply
`kernel_projection.npz` with these native-grid arrays:

```
frequency_hz
quadrature_amplitude_projection_x
quadrature_amplitude_projection_y
c_psd_hz_inv_per_unit_variance_x
c_psd_hz_inv_per_unit_variance_y
```

Projection amplitudes are signed and must have absolute value <= 1.
Their channel assignment, phase and physical normalization require
measurement; no 50:50 or equal X/Y assumption is made. The native C PSD
kernel is an essential measured input, not synthesized in this helper.

```bash
python -m configuration_A.derive_qwp_calibration \
  --kernel-projection-npz /work/kernel_projection.npz \
  --optics-json /work/optics.json \
  --output-calibration-npz /work/ape_calibration.npz
```

Then run `python -m configuration_A.build_response` with the resulting
calibration and explicit provenance, followed by `python -m configuration_A.run`.
The A runner rejects a response file generated by another configuration.
Only APE response preparation is new; vendor signal likelihood remains
`peak_shape=[1]`.
