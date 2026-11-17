# Step 5: APE-specific deterministic axion response for native C-LPSD PSD

**Implementation status:** the source-preserving model interface and its
synthetic tests are complete. The final *measured* optical response, the
physical lab-frame halo, and the native C response calibration are NOT
established by this code. **No physical limits are authorized.**

## Why Step 5 needs APE physics but no vendor-code changes

LIGO `get_beta()` uses strain transfer functions and scalar-DM coupling.
APE detects an axion-induced polarization response. The existing APE
development source already provides a **provisional** Case-2 QWP formula:

```
R(f) = tau_s * unit_conversion * sqrt(
         2*rho_dm_gev4 / ((pi/finesse_s)^2 + 4*sin(pi*f*tau1_s)^2))
```

This is an *amplitude per coupling* (peak-angle convention), not a final
detector transfer measurement. The two APE physics files below are copied
**with identical Git blobs** from `APE-analysis-code`,
`development/likelihood` commit
`3a8f2e94f6b1509fe618b0d25d07e5f53e2dff11`:

- `physics/qwp_signal.py`: `7dcee0f86c8bc314cd46c525b37724ee7c8cb189`
- `physics/axion_halo.py`: `d6273efe62e40e9ee65876c012ab552011824c97`

The optional `maxwell_halo_bin_weights()` is an **unboosted benchmark**,
not the laboratory-frame SHM and not a production coupling template.
Supply an independently validated halo frequency distribution for inference.

## Critical correction: C-D = PSD, C-S = PS

Original `LPSD/IO.c` explicitly maps output column `D` to
`data->psd` and `S` to `data->ps`. APE's historical
`FORMAT "fDNbri"` uses **D**, hence one-sided PSD `rad^2/Hz`.
The earlier development QWP code produces **PS** `rad^2` from a
unit-peak power-window kernel. Substituting that output directly for
`peak_norm` in the **PSD** likelihood is dimensionally wrong.

The new `ape_axion_signal_psd.py` reuses the **unchanged** QWP amplitude
function and requires a supplied **native C PSD response kernel**
`K_D[m,j,k]` in `Hz^-1`, calibrated to the C `D` output from a
unit-time-variance sinusoid in the same window/method configuration.
There is **no invented ENBW or PS-to-PSD multiplier** in the code.

For explicitly supplied halo weights `h_k` (sum 1) and quadrature power
response fractions `a_m`, the deterministic expected signal PSD per
`mu=(g/g_ref)^2` is

```
T_D[m,j] = (g_ref^2 / 2) * a_m *
           sum_k (h[k] * R(f[k])^2 * K_D[m,j,k])
PSD_model[m,j] = B_D[m,j] + mu * T_D[m,j]
```

The `1/2` is the established provisional **peak-to-variance** convention,
not a new correction. `K_D` carries the dimensional PSD normalization.
Neither an extra `sqrt(f/Q_a)` nor a duplicate halo/shape factor is applied.

## Minimal new APE-specific module

`ape_axion_signal_psd.py` provides:

- `build_response_psd(...)`: returns `[2, J]` response on the **exact**
  native C frequency grid, nonzero only in selected output bins. It requires
  explicit halo support/weights, **PSD** kernel, optical parameters,
  X/Y responses, and reference coupling; no defaults are invented.
- `save_provisional_response(...)`: saves exactly the four fields expected
  by the unchanged Step-3 `load_peak_norm()`:
  `frequency_hz`, `peak_norm_x`, `peak_norm_y`, `metadata_json`.
  It marks `PROVISIONAL_NOT_VALIDATED_FOR_PHYSICAL_LIMITS` and refuses
  overwrites.
- `inject_expected_psd(B, T, mu=...)`: deterministic **power-domain**
  synthetic injection; no randomly drawn axion phase/amplitude, no claim that
  it exercises original C-LPSD or proves coverage.

The current Step-4 optimizer bridge uses only `peak_shape=[1]`. Here,
the complete supplied halo and estimator response is **already folded into**
`peak_norm`; do **not** multiply a second halo shape. This allows one-bin
software checks, **not** a correlated multi-bin likelihood or verified
95%-confidence upper limits. Restrict the Step-4 analysis to the response's
supported native-frequency interval.

## Synthetic tests run

From `APE_equivalent_external/`:

```bash
PYTHONPATH=. python -m unittest discover -s ape_post_lpsd/tests \
  -p 'test_ape_axion_signal_psd.py' -v
```

Eight tests pass locally, covering parity with the existing provisional QWP
PS formula given an **explicit** PS-to-PSD response calibration, squared-
coupling scaling, independent X/Y outputs, native grid, invalid halo/kernel
inputs, NPZ interface, deterministic injection, and unboosted halo rejection.
The numerical optical parameters in the tests are **synthetic fixtures**.

## Mandatory remaining physical gates

1. Obtain/derive the correct **measured** cavity/two-QWP amplitude *and phase*
   response, including signal readout projection and peak/RMS conventions.
2. Obtain a physically sourced **laboratory-frame axion halo** and verify
   discretization/coverage. The optional unboosted helper is not sufficient.
3. Measure or source **native C method-2 `D` PSD kernel** for the actual
   APE flags, frequency grid, Kaiser window, overlap and averaging. Check
   its normalization with direct synthetic C tone injections.
4. Determine X/Y noise covariance and **correlations between overlapping
   C-LPSD frequency bins**. Do not assume independent local observations
   for an unvalidated multi-bin likelihood.
5. Validate genuine `findiff` Fisher uncertainty, injection recovery,
   physical coupling conversion and correct coverage on representative data.

Until all relevant gates pass, Step 5 is a **provisional, dimensionally
explicit scientific interface**, not a calibrated axion search result.

## Fidelity

- 68/68 original upstream vendor files remain unchanged.
- 50/50 prior APE preprocessing files remain unchanged.
- All original C, background, likelihood, optimizer, UL and candidate code
  remain unchanged.
- The only new algorithmic glue is the APE-specific PSD response
  packaging and deterministic synthetic injection.
