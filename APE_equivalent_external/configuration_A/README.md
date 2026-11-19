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

## Viper job order (4.2875-hour already calculated C-LPSD)

From `APE_equivalent_external/` on Viper, use the existing branch and run
`git pull`. The shared input spectra are the verified November 17 handoff
at `/ptmp/qarokn/APE_4p28h_PSD_OriginalVendor_PostLPSD_20261117/handoff/`.
The background job reuses these X/Y HDF5 files and creates the summed Q PSD.
The original vendor statistics are unchanged.

1. `sbatch --export=ALL ape_q_full_background_viper.sbatch` (ONLY if no equivalent background job is already running).
2. After `background_Q.json` is complete, submit the 20-bin **5 kHz pilot**:
   `sbatch --export=ALL,APE_LKL_MODE=pilot ape_q_likelihood_viper.sbatch`.
   Inspect `ape_Q_likelihood_*.out`, `*.err` and the pilot's fitted \mu/\sigma
   before scaling up. Scan bounds default to `log10(mu)=-8..8` for
   `mu=(g/g_ref)^2`; never reuse scalar-strain bounds `-45..-32`.
3. After approving the pilot's runtime and fitted values:
   `bash ape_q_submit_full_viper.sh`. It computes the exact number of
   native-grid tasks and submits bounded Slurm array waves.
4. **Only after ALL array tasks complete:** `sbatch --export=ALL
   ape_q_finalize_viper.sbatch`. This verifies full native-grid coverage,
   calls the unchanged original MaxLKL combiner, original expected-UL
   `get_upper_limits_approx` (with the vendor's 1000-point interpolation
   strategy), and original candidate screening/cluster functions.

No observed-limit code is added. The 65 Hz–8 kHz scan includes approximately
4.8 million native frequency bins. Since each bin calls the original
multi-seed optimizer, finishing the FULL scan in a single day cannot be
promised without an actual measured pilot runtime and adequate Viper
parallel resources. An unchanged four-hour data duration can be retained
for smaller, explicitly labeled frequency-band pilot results.

The physical `g` upper-limit curve is provisional until the source PSD
calibration (0.0284 W/V versus older 0.020 W/V) and APE optical signal
normalization are certified against the actual processed input provenance.
