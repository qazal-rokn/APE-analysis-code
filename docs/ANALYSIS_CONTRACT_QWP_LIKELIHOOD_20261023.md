# APE QWP likelihood Analysis Contract (23 October 2026)

## REQUESTED
Integrate APE-only QWP finesse response with the original external vendor likelihood interface on existing development/likelihood; keep main, external code, original calibration, C binary and research data unchanged. All scientific parameter values supplied by user. Update README and equations.

## CURRENT BEHAVIOR
The previous adapter restricted peak_norm to a scalar, while the vendor accepts a frequency-dependent and observation-dependent normalization. The prior QWP physics module was only in an external ZIP.

## PLANNED IMPLEMENTATION
APE adapter delegates the original six-argument vendor stats.log_likelihood without changes to vendor statistics. Factor a physics-generated signal-power template into vendor peak_norm and normalized peak_shape. QWP physics module takes all physical inputs explicitly.

## SCIENTIFIC DEFINITIONS AND UNITS
R(f) = tau * C_unit * sqrt(2*rho_DM / ((pi/F_s)^2 + 4*sin(pi*f*tau1)^2)), nominal high-finesse resonant Case-2 signal amplitude per GeV^-1.

T[m,j] = (g_ref^2/2) * a[m] * sum_k h[k]*R(f_k)^2*W[m,j,k], where h is a supplied normalized halo distribution, W is the supplied native C-LPSD PS kernel, and a is user-supplied quadrature power response. T is rad^2 (native PS). The factor 1/2 assumes R means peak sinusoidal amplitude. No extra sqrt(f/Q_a) term in R. mu=(g/g_ref)^2 for physical nonnegative APE signal. Generic vendor likelihood continues accepting its own signed mu convention.

Vendor Y and bkg are log power (M,N), peak_norm broadcastable to (M,N), peak_shape normalized (N,), model_args (M,3). Do not mix PS rad^2 and PSD rad^2/Hz.

## RISKS AND CONTROLS
Peak/RMS convention, cavity mode/detuning, X/Y response, realistic halo, exact Constant-Q kernel and inference coverage remain physically unvalidated. Overlapping C-LPSD bins are correlated. The helper's unit-sum shape is a vendor interface, not a physical conservation law across bins.

## WILL NOT DO
No new branch, no original vendor edits, no C changes, no real-data reads, no old-result overwrite, no modified cleaning/calibration, no invented physical defaults.

## FILES THAT WOULD CHANGE
APE adapter, APE-only physics file, focused tests, real-vendor checking script's call signature, README, this contract, and EQUATION_REFERENCE.tex.

## VALIDATION PLAN
Only minimal synthetic tests now: original loss-vs-finesse algebra; nonuniform M-by-N peak_norm delegation; power template convention; invalid input. Full true-vendor, C kernel, injection and coverage tests required later before physical inference.

## UNRESOLVED DECISIONS
Actual measured polarization-specific finesse, QWP peak/RMS interpretation, X/Y signal response, discretized halo distribution, measured C-LPSD kernel, and real-data approval.

## REAL-DATA GATE
Not approved. These functions are an integration interface, not a validated coupling-limit pipeline.

## STATUS
User approved narrow development-only work on October 23, 2026.
