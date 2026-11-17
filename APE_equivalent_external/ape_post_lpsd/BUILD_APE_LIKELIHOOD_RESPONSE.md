# APE g-a-gamma-gamma to original vendor PSD likelihood: minimal bridge

This change reuses the EXISTING provisional APE QWP Case-2 amplitude
(`ape_post_lpsd.physics.qwp_signal`) and EXISTING PSD builder
(`ape_post_lpsd.ape_axion_signal_psd.save_provisional_response`).
It adds only a command-line input bridge, without changing numerical methods,
original vendor software, measured PSDs, or background fits.

## Physics convention (unchanged)

`mu = (g_a_gamma_gamma / g_ref)^2`

`D_model[q,j] = B[q,j] + mu * T_D[q,j]`

`T_D[q,j] = (g_ref^2 / 2) * a[q] * sum_k h[k] * R_QWP(f_k)^2 * K_D[q,j,k]`

Here `T_D` is in rad^2/Hz per unit mu and `K_D` must be the independently
validated original C **D-column PSD**, not S-column power spectrum.
The factor 1/2 is the existing peak-to-variance convention, not a new factor.

## Required existing inputs

- Full native-grid `APE_0_15435_X.h5` containing `frequency`.
- Existing `validate_native_c_kernel.export_kernel` format NPZ with keys
  `native_frequency_hz`, `output_bin_indices`, `halo_frequency_hz`,
  `native_c_psd_kernel_per_variance_hz`, `metadata_json`.
- Explicit `physics.json` with exactly these fields:
  `normalized_halo_weights` (correspond to kernel modes and sum to 1),
  `quadrature_power_fraction` (two numbers X and Y),
  `reference_coupling_gev_inv`, `finesse_s`, `tau_s`, `tau1_s`,
  `rho_dm_gev4`, `unit_conversion`, `halo_source`,
  `calibration_version`.

No calibration number, observed halo, or measured window kernel is
silently inferred or supplied. Actual values and sources must be
scientifically checked by the experiment.

## Command

From the `APE_equivalent_external/` package directory:

```bash
python -m ape_post_lpsd.build_ape_likelihood_response \
  --spectrum-hdf5 /path/to/handoff/APE_0_15435_X.h5 \
  --kernel-npz /path/to/verified_native_c_D_kernel.npz \
  --physics-json /path/to/explicit_ape_qwp_parameters.json \
  --output-npz /new/path/ape_gagg_psd_response.npz
```

Output is EXACTLY the existing four-array NPZ understood by
`ape_signal_response.load_peak_norm`: `frequency_hz`,
`peak_norm_x`, `peak_norm_y`, and `metadata_json`.
The original Step-4 likelihood adapter is unchanged and consumes this NPZ.

**Full-grid requirement:** all native output bins must be covered by the
kernel. The `--allow-partial-grid` switch exists only for diagnostics.
No extrapolation, rebinning, or zero-filling of absent kernel coverage is
permitted in the full-band mode.

**Still required before running the actual 4.28-hour likelihood:** a real,
verified native C PSD kernel for the relevant frequency/halo modes and
explicitly reviewed APE QWP/halo parameters. The bridge does not create
such measurements or claim the physics was calibrated. Full-band 5-million
bin x multiple halo-mode kernels may have substantial memory/storage costs.

**Scientific status:** The generated response is explicitly labeled
`PROVISIONAL_NOT_VALIDATED_FOR_PHYSICAL_LIMITS`; the unchanged original
vendor statistics, optical phase, halo, correlation and coverage questions
remain open.

## Focused regression tests

```bash
python -m unittest discover -s ape_post_lpsd/tests \
  -p 'test_build_ape_likelihood_response.py' -v
```

The six tests check bridging, completeness, frequency equality, PSD units,
required physical inputs, and refusal to overwrite. Tests use *synthetic*
data and do not establish empirical optical or C-kernel normalization.

## Source modification scope

- Added `ape_post_lpsd/build_ape_likelihood_response.py`.
- Added `ape_post_lpsd/tests/test_build_ape_likelihood_response.py`.
- Added this documentation.
- **No original vendor or C source was changed.**
- **No old data, background JSON, or compiled artifact was changed.**
