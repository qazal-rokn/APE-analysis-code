# Analysis Contract: Reduced-Block APE Cleaning

## REQUESTED
Use the successful earlier APE cavity-lock cleaning architecture, with one
researcher-approved modification: eta0 quality is no longer a fixed
eta0 >= 0.019 cut. It uses the same local relative-reference philosophy as the
cavity stage.

## APPROVED METHOD

### Reduced records
Raw synchronized data are reduced into non-overlapping 0.05 s records. No
science channel is smoothed or interpolated. Reduced statistics are used only
to make quality decisions.

### Cavity quality
For every reduced record:
- compute cavity-transmission average;
- compute cavity-transmission minimum;
- use the physical midpoint timestamp.

For each record, construct a centered 60 s local reference from eligible
cavity averages using the 90th percentile with linear interpolation at
(n - 1) * q.

A record is provisionally cavity-good when

    cavity_minimum >= 0.90 * local_q90(cavity_average).

Recompute the reference from provisionally good records three times, exactly
following the earlier validated cavity-quality architecture. The final
cavity-good decision uses the last refined reference.

### Preliminary eta0 quality quantity
For each reduced record, compute averages of X2f, Y2f, and cavity transmission.
Use the separately approved eta0 equation

    R2f = hypot(<X2f>, <Y2f>)
    I0  = C_trans * <V_trans>
    eta0 = sqrt(2*sqrt(2)*R2f/(H2*q*G*I0)).

Eta0 records from cavity-bad blocks are not eligible to define the eta0
reference.

### Eta0 quality
For each cavity-good reduced record, construct a centered 60 s local reference
from eligible eta0 values using the same 90th percentile definition.

A record is provisionally eta0-good when

    eta0 >= 0.90 * local_q90(eta0).

Recompute the eta0 reference from provisionally eta0-good AND cavity-good
records three times. The final eta0-good decision uses the last refined
reference.

The final mask is

    final_good = cavity_good AND eta0_good.

### Mapping back to raw samples
Every raw sample in a reduced record inherits that record's quality state.
Original sample indices are preserved. Bad records are never deleted and
neighboring good records are never concatenated across a bad interval.

## PRESERVED FEATURES FROM THE EARLIER VALIDATED CODE
- centered physical-time windows;
- exact timestamp ordering;
- local 90th-percentile reference;
- three reference-refinement iterations;
- cavity average for the reference;
- cavity minimum for the lock decision;
- separate cavity and eta0 masks;
- no interpolation;
- no fallback reference;
- no smoothing;
- no gap filling;
- no padding neighboring records around a rejected record;
- slow drift may track the local reference.

## INTENTIONAL CHANGE FROM THE EARLIER CODE
The earlier PEM-quality module used a fixed

    eta0 >= 0.019

acceptance threshold. The new approved rule replaces it with

    eta0 >= 0.90 * local_q90(eta0)

using the same 60 s centered window and three refinement iterations.

## FULL CALIBRATION GATE
Full rho_X/rho_Y calibration is allowed only where final_good is true.

## VALIDATION PLAN
Synthetic tests must verify:
1. slow cavity and eta0 drifts remain accepted;
2. a reduced cavity block with a deep minimum is rejected even if its average
   remains high;
3. a reduced eta0 glitch is rejected;
4. stable eta0 below 0.019 can pass when it is normal relative to its local
   reference;
5. references remain continuous across raw reader chunk boundaries;
6. the final mask gates rho_X/rho_Y calibration.

Real-data cleaning must be reviewed again after those tests pass.

## STATUS
APPROVED BY RESEARCHER AND IMPLEMENTATION AUTHORIZED
