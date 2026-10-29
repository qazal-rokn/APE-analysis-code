# APE spectral units and historical PSD-only boundary (2026-10-29)

## Scope

Normal native C direct-PS output is already power spectrum [rad²], with its
native estimator normalization. Never introduce an additional f/Q_a,
sqrt(f/Q_a), or coherence-time factor in the normal C/vendor pipeline.
Keep X and Y separate as input observations in the original power-based
vendor background model and likelihood.

## Historical PSD-only files

An archived PSD [rad²/Hz] can be converted to an estimator band power [rad²]
**only** with a separately identified, documented *effective bandwidth*:

    PS [rad²] = PSD [rad²/Hz] * B_eff [Hz].

Likewise an archived ASD [rad/sqrt(Hz)] would require squaring before
multiplying by B_eff. Do not assume column 2 of a C seven-column fDSNbri
output is the only available spectrum: its third column is the native PS.

The value B_eff is NOT automatically f/Q_a. Coherence time and C estimator
effective bandwidth are different questions; using f/Q_a is not authorized
without a proven match to the original estimator's normalization. Native
C PS is preferred when available. This adapter is for a separately named
historical product and must never overwrite or relabel its PSD source.
The historical C_trans=0.020 W/V product must not be labeled as the
newly calibrated C_trans=0.0284 W/V data.

## Amplitude and coupling

The direct amplitude diagnostic is A_XY=sqrt(PS_X+PS_Y) [rad].
Do not pass that amplitude to the original power-domain vendor likelihood,
which treats each quadrature's noise distribution separately.

For physical coupling conversion use *fitted signal* power P_sig [rad²],
not raw observed signal-plus-background power. The corresponding amplitude
is A_sig=sqrt(P_sig) [rad]. With a validated amplitude response
R(f) [rad/(GeV^-1)], g=A_sig/R(f) [GeV^-1].

The existing original vendor parameterization mu=(g/g_ref)^2 is
equivalent: g=g_ref*sqrt(mu) for mu>=0. A multi-quadrature response
requires approved fractions and RMS/peak convention; their physical
validation remains open. A 95% expected UL is not an observed limit.

## Implementation and non-authorization

`src/ape_c_lpsd/analysis/spectral_units.py` provides pure, guarded
unit-boundary functions, without changing C, vendor statistics, existing
pipeline execution, calibration, or the physical model. It DOES NOT
infer unknown bandwidths, convert or load VIPER files, or authorize
real-data likelihood processing. Run its unit tests on the checkout.

Before running the archived-data downstream diagnostic, record: exact
archived paths, schema and C format, PSD versus ASD versus native PS,
effective bandwidth provenance (if PSD-only), X/Y alignment, historical
calibration provenance, and explicitly authorized diagnostic scope.
