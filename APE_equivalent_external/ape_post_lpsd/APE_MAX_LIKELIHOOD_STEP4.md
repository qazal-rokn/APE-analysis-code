# Step 4: Minimal APE input bridge to the ORIGINAL maximum-likelihood optimizer

## Precisely what changed

A new **APE-only module** `ape_max_likelihood.py` calls the *original*
external `max_lkl_executable.calc_max_lkl()`, `fit_background.bkg_model()`,
`LPSDOutput`, `LPSDJSONIO` and `MaxLKLIO` without editing any vendor code.

The original `max_lkl_executable.main()/make_args()` are **not** directly
usable for APE because they construct `peak_norm` from LIGO-specific
`get_A_star()` and `stats.get_beta()`. Rather than patching the external
source, we prepare the exact original optimizer arguments and pass the
pre-existing, explicitly supplied APE `peak_norm` from Step 3.

The fitted statistical model is unchanged:

`PSD_model = exp(bkg) + mu * peak_norm * peak_shape`.

## Intentionally narrow first implementation

- **One quadrature per run**: X or Y. Do *not* assume the two lock-in
  quadratures are independent LIGO interferometers; no joint likelihood.
- **One frequency bin per template**, `peak_shape=[1]`. This is a source-
  parity and single-bin optimizer bridge, **not a verified axion linewidth
  model**. Multi-bin axion shape awaits the separate physical response study.
- Input is one Step-1 `APE_<t0>_<t1>_X.h5` or `_Y.h5`, the Step-2
  original-vendor JSON background fit, and the Step-3
  `frequency_hz/peak_norm_x/peak_norm_y/metadata_json` response NPZ.
- Frequency, log PSD, spline background, skew-normal parameters and response
  values are forwarded with **no numerical changes**.
- The original optimizer performs its signed-`mu` scan, Nelder-Mead search,
  likelihood ratio and Fisher uncertainty; we do not replace these routines.
- Output is unchanged `MaxLKLIO` five-column NPY
  `[frequency, zero_lkl, max_lkl, mu, sigma]`, plus a separately labelled
  provenance JSON. Never overwrite existing results.
- Mu seed-range exponents must be supplied explicitly; LIGO's default
  coupling scale is **not appropriate to assume for APE**. This does not
  redefine the physical parameter.

## Run a bounded validation with your existing spectrum and response

Run from `APE_equivalent_external/`, with all upstream Python dependencies
installed (including `findiff`):

```bash
python -m ape_post_lpsd.ape_max_likelihood \
  --spectrum /real/path/APE_T0_T1_X.h5 \
  --background-json /real/path/background.json \
  --response-npz /real/path/approved_ape_response.npz \
  --output-npy /new/output/ape_X_maxlkl.npy \
  --channel X \
  --ana-fmin 65 --ana-fmax 8000 \
  --min-log10mu YOUR_APPROVED_LOWER_EXPONENT \
  --max-log10mu YOUR_APPROVED_UPPER_EXPONENT \
  --max-frequencies 1
```

Replace the placeholder bounds with numbers appropriate to the **defined
dimensionless** `mu` parameter; do not reuse LIGO's `-40` to `-32`
bounds without an independent scale check. Repeat for Y in a different output
file. Remove `--max-frequencies 1` only after the small-scale tests pass.

The original `max_lkl_executable.main()` still uses LIGO calibration and
must **not** be invoked as the APE analysis command.

## What is and isn't validated

Synthetic local tests (six cases) passed, including **numerical equality with
a direct call to the original `calc_max_lkl()`**, and a guard ensuring LIGO
`stats.get_beta()` is never called. X/Y are processed separately. The exact
original `max_lkl_executable.py`, `stats.py`, `fit_background.py`, and
`LPSDIO.py` source files were used in the local test package.

**Test limitation:** the local environment did not have the vendor's
`findiff` dependency and could not install it. A **test-only, five-point
central second-derivative stand-in** was used to let the original
`sigma_at_point()` exercise its interface; this does **not** certify the
numerical Fisher uncertainty. The committed tests require the genuine
`findiff` package and do not include the substitute.

**Not established:** calibrated APE cavity/QWP axion response, definition
of `mu` in terms of `g_{a\gamma\gamma}`, axion multi-bin spectral shape,
X/Y covariance, absolute coupling limits, candidate significance/coverage,
or parity on real APE spectra. Step 3's response loader checks a user-supplied
normalization but **does not itself derive or certify that physics**.

## Re-run test suite

From `APE_equivalent_external/`, in a full external-program environment:

```bash
python -m unittest discover -s ape_post_lpsd/tests \
  -p 'test_ape_max_likelihood.py' -v
```

## Change budget

- External C/Python vendor source files changed: **0**
- Existing APE pre-LPSD imported files changed: **0**
- Original statistical routines copied or reimplemented: **0**
- New APE-specific code: **one opt-in single-quadrature orchestration file**,
  plus its tests and this explanation.

This is only the minimal **optimizer-input integration**. We must establish
the physical response and multi-bin model, then revalidate the statistical
assumptions before interpreting outputs as axion constraints.
