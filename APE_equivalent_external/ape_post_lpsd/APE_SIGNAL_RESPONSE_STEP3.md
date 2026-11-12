# Step 3: minimal APE response boundary for unchanged vendor likelihood

## Source-level finding

The original external `scalar_dark_matter_lpsd/lvk/max_lkl_executable.py`
computes `peak_norm` in its `make_args` function using
`stats.get_beta(freq_Hz, f_A_star[ifo])`, which requires H1/L1 strain transfer
functions. The original `stats.log_likelihood(mu, Y, bkg, peak_norm, peak_shape,
model_args)` has **no H1/L1 dependency**. It forms the additive one-sided
power spectrum model

`exp(bkg) + peak_norm * peak_shape * mu`.

Therefore, **we must replace the LIGO response at the input boundary, not the
statistical likelihood formula**. The original `max_lkl_executable.main()` is
not yet APE-ready: it still calls `get_A_star()` and `stats.get_beta()`, and
its optimization ranges are tuned to the original signal parameter. This
step does **not** claim otherwise.

## Implemented, and only what is necessary now

`ape_signal_response.py` provides:

1. `load_peak_norm(path, native_frequency_hz)`: reads a separately calibrated
   X/Y signal-power response, checks that its frequencies **exactly match**
   the original C-LPSD frequency grid, and supplies `peak_norm[2, N]`.
2. `vendor_likelihood_inputs(...)`: checks the original vendor argument shapes,
   positive skew-normal scales and normalized `peak_shape`, then returns the
   **same arrays without modifying any numbers**.

It does not resample, add X and Y, modify PSD units, calculate a new
background, alter the external optimizer, or insert a guessed axion response.

## Required input contract

A `.npz` file with four items:

- `frequency_hz` = native C-LPSD frequencies;
- `peak_norm_x`, `peak_norm_y` = nonnegative, one-sided **additive PSD per
  unit of the parameter `mu`** for the measured rho_X and rho_Y channels;
- `metadata_json` = scalar JSON string with `response_units`, `mu_definition`,
  `response_source`, `calibration_version`, and `spectral_shape_convention`.

Use `response_units: "rad^2/Hz per unit mu"`. A metadata string documents the
source but **does not constitute scientific approval**. The response arrays
and their physics must be independently checked before they can be used for
APE sensitivity/upper limits.

**IMPORTANT**: the unchanged vendor likelihood multiplies the response by
`peak_shape`, which must sum to 1. Do not already fold that same peak shape
into `peak_norm`, or the shape will be applied twice. The normalization must
be dimensionally consistent with PSD and the chosen parameter `mu`.

## Physical aspects **not yet implemented or validated**

- The axion-induced polarization rotation per `g_aγγ`, including cavity/QWP
  enhancement, transfer function, and phase conventions;
- The chosen coupling parameter (`mu = g_aγγ^2` or a well-defined normalized
  alternative) and its reference units;
- The axion spectral line shape, its frequency dependence, and the effect of
  bin widths and averaging on the signal PSD;
- Whether the X and Y spectral residuals may be treated as statistically
  independent in the *joint* likelihood, or require additional covariance work;
- Numerical validation with an APE-specific simulated signal injection.

In particular, the vendor `stats.get_beta` remains unchanged and **must not**
be used as the APE response. We have not modified the existing LIGO-specific
`max_lkl_executable.py` or wired it to this response interface yet. That
requires a separate, narrowly scoped Step 4 after the physical response and
normalization have been confirmed.

## Tests

Run from `APE_equivalent_external/ape_post_lpsd/`:

```bash
python -m unittest discover -s tests -p 'test_ape_signal_response.py' -v
```

The tests are synthetic shape, units, validation and no-transformation checks;
**they do not establish the APE coupling normalization**. Scientific model
verification is a separate gate.

## Change budget

- Original external scientific source lines changed: **0**.
- APE pre-LPSD copied files changed: **0**.
- New code: one isolated APE response-input module and its tests/documentation.
