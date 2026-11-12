# Approved APE calibration equations

Status: researcher approved.

The demodulated quadratures X1f, Y1f, X2f, and Y2f are treated as RMS
voltages. The cavity-transmission monitor is converted to optical power using
the run-specific coefficient C_trans.

For each sample:

[
I_0 = C_{\rm trans} V_{\rm trans}
]

[
R_{2f} = \sqrt{X_{2f}^2 + Y_{2f}^2}
]

[
\eta_0 =
\sqrt{
\frac{2\sqrt{2} R_{2f}}
{H_2 q G I_0}
}
]

[
\rho_X =
\frac{X_{1f}}
{\sqrt{2} q G I_0 \eta_0}
]

[
\rho_Y =
\frac{Y_{1f}}
{\sqrt{2} q G I_0 \eta_0}
]

H2 is the measured photodiode/electronics response at 2f. X and Y science
quadratures remain separate. The calibration stage performs no smoothing,
interpolation, resampling, phase correction, gain-imbalance correction,
leakage correction, or X/Y combination.

A sample is invalid if the required raw inputs are non-finite, I0 is not
strictly positive, R2f is not strictly positive, or the resulting calibration
denominator is not finite and strictly positive. Invalid samples retain their
original indices and calibrated outputs are NaN.
