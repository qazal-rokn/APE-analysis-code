# Approved APE cleaning definitions

Status: researcher approved.

Cleaning is staged and one-sided. Its purpose is to reject sharp downward
failures while preserving slow drift.

The cavity-transmission reference is the centered 60 s local median of 0.05 s
block medians. A raw sample is cavity-good when

\[
V_{\rm trans}(t) \ge 0.90\,V_{\rm med}(t).
\]

Only cavity-good samples are used to form preliminary

\[
I_0=C_{\rm trans}V_{\rm trans},
\qquad
\eta_0=
\sqrt{
\frac{2\sqrt{2}\,R_{2f}}
{H_2 q G I_0}
}.
\]

The eta0 reference is the centered 60 s local median of 0.05 s eta0 block
medians, using cavity-good samples only. A sample is eta0-good when

\[
\eta_0(t) \ge 0.90\,\eta_{0,\rm med}(t).
\]

The final mask is

\[
M_{\rm final}=M_{\rm cavity}\land M_{\eta}.
\]

There is no upper-side cut. No sample is deleted or concatenated across a bad
interval. Slow drift is intentionally retained. Full rho_X/rho_Y calibration
is performed only where the final mask is true.
