# Analysis Contract: APE Cavity and Eta0 Cleaning

## STATUS

Researcher approved. Implementation authorized.

## CAVITY STAGE

Raw synchronized data are partitioned into non-overlapping 0.05 s records.
For each record (k),

    cavity_average[k] = mean(V_trans)
    cavity_minimum[k] = min(V_trans).

The cavity reference is the centered 60 s 90th percentile of eligible cavity
block averages, with linear interpolation at `(n - 1) * q`.

A block is cavity-good when

    cavity_minimum[k] >= 0.90 * cavity_reference[k].

The reference is recomputed from provisionally cavity-good blocks three times,
exactly following the earlier validated APE cavity-quality architecture.

The cavity stage is unchanged by the eta0 optimization.

## SAMPLE-LEVEL ETA0

Inside a cavity-good candidate segment, eta0 is calculated at the original raw
sample rate:

    I0[n] = C_trans * V_trans[n]
    R2f[n] = hypot(X2f[n], Y2f[n])
    eta0[n] = sqrt(2*sqrt(2)*R2f[n] / (H2*q*G*I0[n])).

For each 0.05 s block,

    eta0_mean[k] = mean(eta0[n] in block k)
    eta0_min[k]  = min(eta0[n] in block k).

The mean is the reference statistic. The minimum is the quality-decision
statistic.

## PROTECTED ETA0 REFERENCE

Within each cavity-good segment, the eta0 reference is past-only:

    eta_ref[k] =
        Q0.90({
            eta0_mean[j] :
            t_k - 60 s <= t_j < t_k,
            eta0_good[j]
        }).

The current block never contributes to its own reference.

The decision is

    eta0_min[k] >= 0.80 * eta_ref[k].

Rejected eta0 blocks do not update the future reference. If a sustained eta0
failure leaves no accepted block inside the preceding 60 s, the last valid
healthy reference is held until eta0 recovers.

## CAVITY-SEGMENT RESET RULE

Eta0 reference state is intentionally reset at the beginning of every
cavity-good segment. The first valid eta0 block mean in that segment bootstraps
the protected reference.

This reset rule is an explicit researcher-approved methodological choice. File
boundaries alone do not reset eta0 state: a cavity-good segment may span
multiple declared-continuous DXD files.

## EFFICIENT LONGEST-SEGMENT SEARCH

The production selection workflow is:

1. read only cavity transmission across the full acquisition;
2. construct the full cavity mask;
3. identify and rank cavity-good segments by duration;
4. starting with the longest candidate, read the full five raw channels only
   over that candidate's sample range;
5. calculate sample-level eta0 and the protected eta0 mask inside that
   candidate;
6. retain the longest surviving eta0-good subsegment;
7. continue to the next cavity candidate only while its full cavity duration
   can still exceed the current best surviving duration.

The search stops once

    current_best_final_duration >= next_unprocessed_cavity_duration.

This stopping rule is exact under the approved cavity-segment reset convention
because eta0 cleaning can only shorten a cavity-good candidate.

## FINAL MASK AND DATA HANDLING

Within an evaluated cavity-good candidate,

    final_good = cavity_good AND eta0_good.

Every raw sample inherits its 0.05 s block quality state. Bad samples are not
deleted, interpolated, smoothed, padded, filled, or concatenated across bad
intervals.

Full rho_X/rho_Y calibration is allowed only where final_good is true.

## VALIDATION REQUIREMENTS

Synthetic validation must verify:

1. slow cavity and eta0 drift remain accepted;
2. cavity block minimum rejects a sharp cavity glitch;
3. sample-level eta0 block minimum rejects a sharp eta0 glitch;
4. a sustained eta0 collapse cannot redefine the protected reference;
5. eta0 recovery is accepted;
6. stable eta0 below the historical fixed 0.019 value can pass;
7. cavity-only planning reproduces the cavity stage of exhaustive cleaning;
8. eta0 candidate cleaning resets and bootstraps at a cavity-good segment start;
9. the exhaustive segmented eta0 result agrees with independent candidate
   cleaning for the same cavity-good segment;
10. final_good gates rho_X/rho_Y calibration.

The earlier 42 h diagnostic with the protected 80% eta0 cut found approximately
288.40 s of eta0-only rejection. File 0012 remained fully accepted and the
known long eta0 failure in file 0003 contributed approximately 286 s.
