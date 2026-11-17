# APE post-LPSD software pipeline: original algorithms, minimal detector-specific interfaces

**Status:** Software connection for separate X/Y single-bin analysis, candidate *screening* and **median expected** asymptotic limits implemented. **Not a validated APE physical detection/exclusion pipeline.** The original vendor source and existing APE pre-LPSD sources are unchanged.

## Recovered existing APE physics

The archived APE discussion and the existing `development/likelihood` branch contain the **provisional** cavity/two-QWP amplitude model:

`R(f) = tau_s * C_unit * sqrt(2*rho_DM / ((pi/F_s)^2 + 4*sin(pi*f*tau1_s)^2))`.

The model and the optional unboosted Maxwellian benchmark are **already copied byte-for-byte** under `ape_post_lpsd/physics/qwp_signal.py` and `ape_post_lpsd/physics/axion_halo.py`. Do not invent another optical model. Both the actual APE complex optical response and the native C-LPSD D-column signal PSD kernel still need measured/verified normalization. Earlier APE lab-halo/Kaiser tests are benchmarks, not full physical certification.

## Original external functions used unchanged

- `LPSDIO.LPSDOutput`, `LPSDJSONIO`, `MaxLKLIO`;
- `lvk.fit_background.fit_background_in_file` and `bkg_model`;
- `lvk.stats.log_likelihood`, `lvk.max_lkl_executable.calc_max_lkl`;
- `lvk.combine_maxlkl_output.main` for chunk concatenation/sorting;
- `lvk.get_candidates.cluster` for cluster selection;
- `lvk.calc_upper_limit.get_upper_limits_approx` for **median expected**, not observed, asymptotic upper limits.

The original vendor `get_candidates.main` cannot be reused directly: it loads LIGO `get_A_star()`, assumes H1/L1 consistency and requires signals in *both independent detectors*; it also contains `pdb.set_trace()`. The original `calc_upper_limit.main` assumes LIGO strain/coupling conversion, comparison datasets and Cardiff filesystem paths. We bypass those specific `main` entry points, leaving the reusable original functions intact.

## New APE-only glue in this change

1. `ape_post_lpsd/ape_max_likelihood.py` adds **only** `--start-index` to the already existing Step 4 APE interface so jobs may scan nonoverlapping native frequency chunks, using the original optimizer unchanged. This is a small patch to *our APE interface*, never the vendor file. It logs the actual native-bin start/end.
2. `ape_post_lpsd/ape_combine_parts.py` verifies channel identity, response/background configuration and chunk-index continuity; then delegates numerical assembly to **unchanged vendor** `combine_maxlkl_output.main()`. It writes a separate provenance JSON and does not overwrite.
3. `ape_post_lpsd/ape_finalize_inference.py` accepts a `MaxLKLIO` NPY and matching APE provenance, reproduces the **original** positive-`mu` q0 screening and Bonferroni convention, delegates clustering to vendor `cluster()` unchanged (with one inactive sentinel to close a cluster at the final bin), and delegates median expected-limit solving to `get_upper_limits_approx()` unchanged. It never calls `get_beta`, `get_A_star` or LIGO-specific two-detector consistency functions.

A physical `g` column is reported *conditionally* if metadata explicitly defines `mu=(g_a_gamma_gamma/g_ref)^2`; it is `g_ref * sqrt(mu_UL)` with **no** LIGO transfer conversion. If `mu` is arbitrary, the g columns are NaN. Provisional status is always clearly recorded.

## Commands from `APE_equivalent_external/`

Install **original vendor Python requirements**, including real `findiff`, in the environment. The local regression tests used a test-only stand-in for its import.

Use your existing pre-LPSD pipeline and Step 1 adapter to create aligned native-grid `APE_T0_T1_X.h5` and `APE_T0_T1_Y.h5`, then the **original** background fitter independently per channel:

```bash
python -m scalar_dark_matter_lpsd.lvk.fit_background --data-path /work/APE_T0_T1_X.h5 --output-path /work/background.json --ana-fmin 65 --ana-fmax 8000
python -m scalar_dark_matter_lpsd.lvk.fit_background --data-path /work/APE_T0_T1_Y.h5 --output-path /work/background.json --ana-fmin 65 --ana-fmax 8000
```

Supply the **already implemented** Step 5 APE response NPZ with correct native frequencies and one-sided PSD units. It is currently provisional and should **not** be used for physical coupling limits until the optical response, C kernel, halo and statistical assumptions are validated.

For X, execute the original optimizer through the existing APE bridge, in consecutive 100-bin example chunks:

```bash
python -m ape_post_lpsd.ape_max_likelihood --spectrum /work/APE_T0_T1_X.h5 \
  --background-json /work/background.json --response-npz /work/ape_response.npz \
  --channel X --output-npy /work/parts/ape_X_0_result.npy \
  --ana-fmin 65 --ana-fmax 8000 --min-log10mu MIN_EXP --max-log10mu MAX_EXP \
  --start-index 0 --max-frequencies 100
python -m ape_post_lpsd.ape_max_likelihood --spectrum /work/APE_T0_T1_X.h5 \
  --background-json /work/background.json --response-npz /work/ape_response.npz \
  --channel X --output-npy /work/parts/ape_X_1_result.npy \
  --ana-fmin 65 --ana-fmax 8000 --min-log10mu MIN_EXP --max-log10mu MAX_EXP \
  --start-index 100 --max-frequencies 100
```

Replace `MIN_EXP`, `MAX_EXP` with scan bounds appropriate to the explicitly defined `mu`; do not use original LIGO defaults blindly. Continue chunk offsets until **all analysis bins** are processed, without overlapping jobs. The two commands above are **only an example**, not an entire scan. Repeat with `--channel Y` and different files.

Then combine using the **original vendor combiner** via the APE provenance guard:

```bash
python -m ape_post_lpsd.ape_combine_parts \
  --results-dir /work/parts --prefix ape_X --output-npy /work/combined_X.npy
```

Note: the combiner checks continuity of supplied chunks but cannot prove that the intended final frequency endpoint was reached; verify the complete requested scan separately.

Screen and calculate **median expected** upper limits in the original vendor statistical convention:

```bash
python -m ape_post_lpsd.ape_finalize_inference \
  --maxlkl-npy /work/combined_X.npy \
  --step4-manifest-json /work/combined_X.json \
  --output-prefix /work/APE_X_diagnostic
```

Files written: `APE_X_diagnostic_screened_candidates.csv`, `APE_X_diagnostic_expected_ul.csv`, and `APE_X_diagnostic_analysis_status.json`. If only one unchunked Step-4 output exists, pass it directly with its own sidecar; do not combine. Repeat all post-LPSD stages independently for Y.

## Validated here / not validated here

**Eleven local synthetic boundary tests pass**. They check delegation to unchanged original functions via mocks, candidate closing at the final bin, expected-limit metadata mapping, original combiner invocation, no accidental X/Y combination, chunk provenance, and non-overwrite behavior. They **do not retest the original vendor numerical methods**. The local environment lacked `findiff`; a test-only stand-in was used for module import and was **not committed**.

**Scientific gates remaining:** real native C-LPSD PSD response-kernel measurement and original-binary provenance; measured cavity/two-QWP response including phase, orientation and peak/RMS convention; physically justified lab-frame axion halo and any full multi-bin correlated likelihood; independent X/Y and overlapping C-bin covariance characterization; APE candidate veto and trials calibration; Fisher uncertainty with genuine `findiff`, injection recovery and APE coverage. The vendor upper-limit helper returns the *median expected limit*, **not a measured observed 95% limit**. Do not relabel it.

**Fidelity:** Original vendor C, background, likelihood, optimizer, clustering and upper-limit function source files remain byte-for-byte unchanged. Existing APE calibrated pre-LPSD source remains unchanged. No original method has been reimplemented to make the software appear complete.
