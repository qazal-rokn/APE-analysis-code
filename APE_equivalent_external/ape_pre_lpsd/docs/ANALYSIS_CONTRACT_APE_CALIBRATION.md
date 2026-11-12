# Analysis Contract: Canonical APE Calibration

## REQUESTED
Build the calibration stage after validated Dewesoft ingestion.

## CURRENT BEHAVIOR
The project reads synchronized raw X1f, Y1f, X2f, Y2f, and cavity-transmission
channels in bounded chunks. A legacy review-only calibration exists separately.

## UNDERSTANDING
The new canonical calibration must use the researcher-approved RMS equations
and must keep rho_X and rho_Y separate.

## PLANNED IMPLEMENTATION
For each raw APEDXDChunk, compute I0, R2f, eta0, rho_X, rho_Y, and a valid mask.
Invalid samples retain their original indices and calibrated outputs are NaN.

## SCIENTIFIC DEFINITIONS AND UNITS
I0 = C_trans * V_trans, in W.

R2f = hypot(X2f, Y2f), in V RMS.

eta0 = sqrt(2*sqrt(2)*R2f/(H2*q*G*I0)), dimensionless.

rho_X = X1f/(sqrt(2)*q*G*I0*eta0), dimensionless.

rho_Y = Y1f/(sqrt(2)*q*G*I0*eta0), dimensionless.

The recorded demodulated quadratures are treated as RMS voltages.
H2 is the measured photodiode/electronics response at 2f.

## RISKS AND CONTROLS
The main risk is an RMS/peak or H2 convention error. The equations are frozen
explicitly and validated with synthetic forward-generated data.

## WILL NOT DO
No smoothing, filtering, interpolation, resampling, clipping, phase correction,
gain-imbalance correction, leakage correction, X/Y combination, cleaning cut,
spectral estimation, or C-LPSD execution.

## FILES THAT WOULD CHANGE
- src/ape_c_lpsd/calibration.py
- src/ape_c_lpsd/ingestion/models.py
- src/ape_c_lpsd/ingestion/__init__.py
- src/physics/APE_CALIBRATION.md
- docs/EQUATION_REFERENCE.tex
- tests/test_calibration.py
- documentation describing the canonical calibration

## VALIDATION PLAN
Forward-generate synthetic raw voltages from known I0, eta0, rho_X, and rho_Y,
then require the calibration to recover the injected values to float64
precision. Also verify invalid physical samples remain indexed and become NaN.

## UNRESOLVED DECISIONS
None for the equations implemented here.

## REAL-DATA GATE
Synthetic validation must pass before applying the new calibration to real APE
data.

## STATUS
APPROVED BY RESEARCHER AND IMPLEMENTATION AUTHORIZED
