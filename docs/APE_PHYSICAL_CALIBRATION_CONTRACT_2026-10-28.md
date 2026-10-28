# APE physical-calibration decisions and outstanding gates (2026-10-28)

This records the researcher's declared calibration values for **future physical APE data processing**. It does not alter the frozen synthetic 5 kHz benchmark, reprocess raw data, or certify real-data axion exclusions. The current existing calibrated HDF5 uses historical values and must remain immutable.

## Confirmed values for the forthcoming real-data calibration contract

| Quantity | Symbol / code parameter | Value | Status |
|---|---|---:|---|
| Transmission photodiode calibration | `C_trans` | **0.0284 W/V** | explicitly reconfirmed by researcher |
| Cavity finesse for the QWP response | `finesse_s` / `mathcal F` | **1500** | explicitly reconfirmed by researcher |
| Photodiode responsivity | `q` | 0.7 A/W | prior approved APE convention |
| Transimpedance gain | `G` | 1e6 V/A | prior approved APE convention |
| QWP response time parameter | `tau_s` | 2.5e-9 s | prior approved Case-2 convention |
| QWP cavity delay parameter | `tau1_s` | 1e-8 s | prior approved Case-2 convention |
| Local halo energy density | `rho_dm_gev4` | 3.1e-42 GeV^4 | prior benchmark convention, approximately 0.4 GeV/cm^3 |
| SI-to-natural time conversion | `unit_conversion` | 1.519e24 GeV^-1 s^-1 | prior convention |
| Axion coherence quality | `Q_a` | 1e6 | prior analysis convention |
| Reference coupling in vendor power model | `reference_coupling_gev_inv` | 1e-10 GeV^-1 | software-normalization reference, not a physical assumption |

**Unconfirmed for final physical inference:** `H_2`, the *measured* transfer/attenuation of the 2f photodiode/electronics at the relevant frequency. Existing HDF5 metadata says `second_harmonic_pd_response = 0.7` but that value alone is not proof of a measured response. The actual X/Y axion-induced quadrature power fractions and phase transfer are also provisional; do not silently assume `alpha_X=alpha_Y=0.5`. The QWP response used in synthetic tests is an amplitude model, not a validated full complex optical response.

## Approved conversion convention (RMS quadratures, separate X/Y)

```text
I0(t) = 0.0284 [W/V] * cavity_transmission_V(t)
R2f(t) = hypot(X2f_RMS(t), Y2f_RMS(t))
eta0(t) = sqrt( 2 sqrt(2) R2f(t) / [H2 q G I0(t)] )
rho_X(t) = X1f_RMS(t) / [sqrt(2) q G I0(t) eta0(t)]
rho_Y(t) = Y1f_RMS(t) / [sqrt(2) q G I0(t) eta0(t)]
```

All parameters and conventions apply at sample level; do not take median eta0 over 60 seconds and do not combine X/Y at this stage.

Existing checked-in physics function `src/physics/qwp_signal.py` expects the following **provisional** QWP amplitude response, with finesse set to 1500 in a *new run-specific* configuration:

```text
R_QWP(f) = tau_s * C_unit * sqrt(2 rho_DM)
         / sqrt( (pi/1500)^2 + 4 sin(pi f tau1_s)^2 )
```

The calibrated native C result is **PS [rad^2]**, not PSD [rad^2/Hz]. A deterministic halo power template already includes native C frequency-window response. In the adapter, `mu = (g/g_ref)^2` and `g_UL = g_ref * sqrt(mu_UL)`, with the vendor's reported limit explicitly **median expected**, not observed.

**Do not add any external `sqrt(f/Q_a)` or `sqrt(Q_a/f)` factor** (or frequency-bin-width/sensitivity multiplier) to the QWP response or coupling conversion; the coherence-matched native C-LPSD transformation handles the relevant integration window. Do not conflate `f` (signal frequency) with `F` (cavity finesse). This is the agreed model convention, not proof that the physical response has been experimentally validated.

## Historical alternative: analytic shot-noise sensitivity forecast (not the current measured-PS likelihood)

Earlier Case-2 sensitivity forecasts included:
- `Ts2 = 2980e-6/4` = 745 ppm, `Tp2 = 2980e-6/4` = 745 ppm,
- `T1s2 = 600e-6` = 600 ppm, `T2s2 = 580e-6` = 580 ppm,
- `P0 = 0.17 W`, `phi_fixed = 90 deg`, `sigma_squared = 1.2e-7`,
- elementary charge `e = 1.60217663e-19 C`.

These values belong to a **separate optical/shot-noise forecasting calculation**. When real `I0(t)`, `eta0(t)`, calibrated X/Y, and measured C spectra are used, do not insert these forecasting inputs as extra multiplicative signal normalization factors.

## Existing immutable calibrated HDF5

`/work/lpsd-project/LPSD/ape_frozen_calibrated_segment.h5`
contains 308717000 synchronized samples at 20000 Hz, with `i0_w` [W], `eta0` [dimensionless], and separate `rho_x`, `rho_y` [rad]. Its stored **C_trans = 0.020 W/V**, `q=0.7`, `G=1e6`, `H2=0.7`; no filtering, interpolation, detrending, resampling or quadrature combination is recorded.

**It is a historical calibration product, not compliant with the newly confirmed 0.0284 W/V contract. Do not edit it in place or pass it as fully approved physical data.**

With fixed raw signal and identical `H2,q,G`, the analytic ratio `rho_new/rho_old = sqrt(0.020/0.0284)` holds. This is **not an instruction to silently rescale PS or g limits**; any corrected real-data product must be independently named, provenance-bearing, and explicitly run-scoped.

## Next bounded step

Create a **metadata-only, fail-closed calibration gate** comparing a new run-specific contract against the HDF5. Require correct `C_trans`, sampling rate, units, exact X/Y separation, preserved provenance, and physical `H2` provenance. Keep this gate separate from synthetic benchmark acceptance. Do not execute C or regenerate raw data as part of the metadata check. Resolve the optical X/Y response before claiming physical coupling limits.
