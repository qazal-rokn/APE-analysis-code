# Configuration A — minimal original-vendor analysis

**One working branch:** `experiment/ape-configuration-a-summed-quadratures-20261119`.

Configuration A has only two Python modules:

- `combine_quadrature_psd.py`: the necessary APE-only handoff that computes
  `PSD_Q = PSD_X + PSD_Y` in **linear** PSD units, then
  `logPSD_Q = log(PSD_Q)`. It preserves the native C grid and produces an
  HDF5 file compatible with the **original** `LPSDOutput`.
- `run.py`: a small CLI that calls the **original**
  `scalar_dark_matter_lpsd.lvk.fit_background` for Q, then passes Q into
  the existing APE-only detector-response adapter
  `ape_post_lpsd.ape_max_likelihood`, which invokes the **original**
  `max_lkl_executable.calc_max_lkl` and `stats.log_likelihood`.

The original vendor `LPSD/` C source, background fitter, likelihood,
optimizer, and upper-limit/statistical modules are unchanged. The model
uses the original `peak_shape=[1]` single-bin convention. No separate
Monte Carlo/halo/window fit is used in Configuration A.

## 1. Combine X/Y and calculate a NEW background using original vendor code

From `APE_equivalent_external/`:

```bash
python -m configuration_A.run background \
  --x-h5 /work/APE_0_15435_X.h5 \
  --y-h5 /work/APE_0_15435_Y.h5 \
  --output-dir /work/combined_Q \
  --ana-fmin 65 --ana-fmax 8000
```

This produces `/work/combined_Q/APE_0_15435_Q.h5` and
`/work/combined_Q/background_Q.json`. Previous X/Y backgrounds are
not valid for Q and are not modified. The vendor fitter uses its
unchanged defaults unless optional `--n-processes` or
`--segment-size` is explicitly supplied. Refit Q before likelihood.

## 2. Generate the APE response using the EXISTING two-QWP formula

The **only** model-specific substitution is `peak_norm_Q = beta_APE_Q`.
The existing `ape_post_lpsd.physics.qwp_signal.qwp_amplitude_per_coupling`
provides the nominal peak response per coupling:

```
R(f) = tau_s * C * sqrt(2*rho_DM / ((pi/finesse_s)^2 + 4*sin(pi*f*tau1_s)^2))
```

With `mu=(g/g_ref)^2`, the APE-specific peak normalization is:

```
beta_APE_Q(f) = 0.5 * [g_ref * R(f)]^2 / (pi * f * 1e-6)
```

The **1e-6 fractional linewidth, pi and single-bin peak_shape=[1]**
follow the ORIGINAL vendor `stats.get_beta` convention, not a new
Monte Carlo/halo template. The 0.5 converts the nominal peak angle
amplitude to mean-square angle. The summed X+Y spectrum is assumed to
recover all calibrated orthogonal angle power; there is no arbitrary
X/Y split.

Previously discussed APE two-QWP working defaults: `finesse_s=1500`,
`tau_s=2.5e-9 s`, `tau1_s=1e-8 s`,
`rho_dm_gev4=3.1e-42 GeV^4`, `unit_conversion=1.519e24`;
`g_ref=1e-10 GeV^-1`. THESE VALUES WERE PROVISIONAL and must not
be passed off as independently confirmed measurements.

The measured APE calibration constants are preserved in response
metadata for provenance only: `I0=0.0284*V_trans W`,
`q=0.7 A/W`, `G=1e6 V/A`, 2f response `0.7`.
They are **not applied again**, because existing native-C inputs
`rho_x` and `rho_y` were already calibrated in radians.

```bash
python -m configuration_A.run response \
  --spectrum /work/combined_Q/APE_0_15435_Q.h5 \
  --output-npz /work/combined_Q/ape_Q_response.npz
```

The command uses those historic defaults. To use final approved optical
values instead, specify `--finesse-s`, `--tau-s`, `--tau1-s`,
`--rho-dm-gev4`, `--unit-conversion`, and/or
`--reference-coupling-gev-inv`. No original science code changes.
The new response NPZ has only `frequency_hz`, `peak_norm_q`,
`metadata_json`, on the exact original native C grid.

## 3. Run original vendor likelihood on the summed Q spectrum

```bash
python -m configuration_A.run likelihood \
  --spectrum /work/combined_Q/APE_0_15435_Q.h5 \
  --background-json /work/combined_Q/background_Q.json \
  --response-npz /work/combined_Q/ape_Q_response.npz \
  --output-npy /work/combined_Q/likelihood_Q.npy \
  --ana-fmin 65 --ana-fmax 8000 \
  --min-log10mu MIN --max-log10mu MAX
```

Select scan bounds appropriate to this now explicitly defined
`mu=(g/g_ref)^2` and use the unchanged original optimizer.

**State:** code updated on GitHub; no real-data Q background refit or
full-band likelihood execution has been performed here. The provisional-response generator now exists, but historical optical values and absolute signal-normalization conventions must still be approved for physical coupling results.
Configuration B stays separate; existing exploratory scripts elsewhere
in the repository are not imported by these two modules.

## 4. Upper limits from the original five-column likelihood output

The existing APE inference adapter imports and reuses the **unchanged
vendor** `get_upper_limits_approx()` for the 95% median **expected** limit.
The upstream program has **no observed-limit function**. To avoid
misrepresenting this, Configuration A now includes an explicitly labeled
APE-only **observed asymptotic** calculation. It uses the original vendor
signed best fit `mu_hat` and Fisher `sigma`, and the same asymptotic
`q_mu` CDF convention:

```
mu_UL,observed = max(0, mu_hat + Phi^{-1}(0.95)*sigma)
g_UL,observed  = g_ref*sqrt(mu_UL,observed)
```

A strongly negative `mu_hat` that gives no positive Wald upper bound
is flagged with `zero_boundary_flag=1`. The observed result is an
**asymptotic approximation, not a complete profile-likelihood inversion or
a function copied from the upstream code**. Its confidence coverage has
not been calibrated for APE.

After the Q likelihood produces a five-column `likelihood_Q.npy` and
a sidecar `likelihood_Q.json`, run:

```bash
python -m configuration_A.run limits \
  --maxlkl-npy /work/combined_Q/likelihood_Q.npy \
  --step4-manifest-json /work/combined_Q/likelihood_Q.json \
  --output-prefix /work/combined_Q/ape_Q
```

Outputs: `ape_Q_expected_ul.csv` (unchanged vendor helper),
`ape_Q_observed_ul_asymptotic.csv` (new APE-only approximation),
`ape_Q_screened_candidates.csv` (vendor-style screening), and
`ape_Q_analysis_status.json` (provenance). No vendor files changed.
