# Analysis Contract: Reduced-Block APE Cleaning

## REQUESTED
Keep the successful earlier APE cavity-lock cleaning architecture unchanged.
For eta0, replace the historical fixed threshold and the earlier centered
relative-reference trial with the researcher-approved protected-reference
method validated on the 42 h APE data set.

## APPROVED METHOD

### Reduced records
Raw synchronized data are partitioned into non-overlapping 0.05 s records.
Reduced statistics are used only for quality decisions. Original sample
indices are preserved and no bad interval is deleted or bridged.

### Cavity quality: unchanged
For every reduced record:
- compute cavity-transmission average;
- compute cavity-transmission minimum;
- use the physical midpoint timestamp.

For each record, construct a centered 60 s local reference from eligible
cavity averages using the 90th percentile with linear interpolation at
`(n - 1) * q`.

A record is provisionally cavity-good when

    cavity_minimum >= 0.90 * local_q90(cavity_average).

Recompute the cavity reference from provisionally good records three times,
following the earlier validated cavity-quality architecture. The final
cavity-good decision uses the last refined reference.

### Sample-level eta0
After the cavity stage is defined, calculate eta0 at the original raw-sample
rate from the recorded second-harmonic quadratures and cavity transmission:

    I0[n]  = C_trans * V_trans[n]
    R2f[n] = hypot(X2f[n], Y2f[n])
    eta0[n] = sqrt(2*sqrt(2)*R2f[n]/(H2*q*G*I0[n])).

For each 0.05 s block k, calculate

    eta0_mean[k] = mean(eta0[n] in block k)
    eta0_min[k]  = min(eta0[n] in block k).

A block with any invalid eta0 sample has invalid eta0 block statistics and
cannot be eta0-good.

### Protected eta0 reference
The eta0 reference is past-only. For block k, use the 90th percentile of
previously accepted eta0 block means from the preceding 60 s:

    eta_ref[k] =
        Q0.90({
            eta0_mean[j] :
            t_k - 60 s <= t_j < t_k,
            cavity_good[j],
            eta0_good[j]
        }).

The current block never contributes to its own reference.

A cavity-good block is eta0-good when

    eta0_min[k] >= 0.80 * eta_ref[k].

Only accepted cavity-good AND eta0-good block means are inserted into the
future eta0 reference. Rejected eta0 blocks therefore cannot drag the baseline
downward and redefine a sustained failure as normal.

If a long failure or cavity-bad interval empties the 60 s accepted-history
window, hold the last valid healthy eta0 reference until a block passes again.

At the beginning of the complete acquisition, before previous healthy history
exists, bootstrap the reference from the first valid cavity-good eta0 block
mean. This matches the validated real-data diagnostic and assumes acquisition
begins in its normal operating state.

For the efficient longest-segment workflow, cavity-good segments are treated as
independent eta0 candidates. The protected eta0 reference is reset at the start
of each cavity-good candidate segment and bootstrapped from its first valid
eta0 block. File boundaries alone do not create a reset when a cavity-good
segment crosses between files.

### Final mask
The final reduced-record mask is

    final_good = cavity_good AND eta0_good.

Every raw sample in a reduced record inherits that record's quality state.

## PRESERVED FEATURES
- cavity cleaner unchanged;
- cavity centered physical-time windows;
- cavity local 90th-percentile reference;
- cavity three-refinement procedure;
- cavity average for reference and minimum for decision;
- eta0 uses the same 0.05 s block duration;
- separate cavity and eta0 masks;
- no interpolation;
- no smoothing;
- no gap filling;
- no deletion;
- no concatenation across rejected intervals;
- one-sided lower quality cuts only;
- slow eta0 drift can update the protected reference while data remain healthy.

## INTENTIONAL ETA0 CHANGES
The old historical PEM-quality module used a fixed

    eta0 >= 0.019

threshold.

The first relative-reference implementation instead used one eta0 value derived
from block-averaged X2f/Y2f/Vtrans and a centered 60 s q90 reference with three
refinement passes. Real data exposed a failure mode: a several-minute eta0
collapse in file 0003 eventually dragged that local reference down and began to
look normal.

The approved replacement is therefore:
1. sample-level eta0 calibration;
2. 0.05 s eta0 mean for the reference statistic;
3. 0.05 s eta0 minimum for the decision statistic;
4. past-only 60 s q90 from accepted history;
5. failed eta0 records never update the reference;
6. hold the last healthy reference when the accepted-history window empties;
7. one-sided lower decision at 80% of the protected reference;
8. reset/bootstrap the protected eta0 state at each cavity-good candidate
   segment for the efficient longest-segment search.

## EFFICIENT LONGEST-SEGMENT SEARCH

For the specific task of finding the single longest final-good continuous
segment, the researcher approved an optimized search that does not calculate
sample-level eta0 over the entire acquisition.

The search is:

1. run the unchanged cavity cleaner across the full declared-continuous
   acquisition;
2. identify and rank cavity-good segments by duration;
3. evaluate eta0 only inside cavity-good candidate segments that can still beat
   the best final-good segment found so far;
4. stop once the next unprocessed cavity-good segment is no longer than the
   current best final-good result.

For this optimized search only, the protected eta0 history is intentionally
reset at the start of each cavity-good candidate segment and bootstrapped from
that segment's first valid eta0 block mean. This is an explicit
researcher-approved methodological choice for segment selection, not a hidden
implementation detail.

The eta0 equations, 0.05 s mean/minimum statistics, past-only q90, 80% lower
cut, rejected-block exclusion, and held-reference behavior within each
candidate segment remain unchanged.

## EFFICIENT LONGEST-SEGMENT SEARCH

For the production longest-segment search, the approved computational order is:

    full 42 h -> cavity-only cleaning -> rank cavity-good segments
              -> eta0 cleaning only on candidates that can still win.

For this efficient search only, the protected eta0 history is intentionally
reset at the start of each cavity-good candidate segment. The first valid eta0
block mean inside that candidate bootstraps the reference. This is an explicit
researcher-approved change from carrying eta0 history continuously across the
entire acquisition.

Candidate segments are processed from longest to shortest. After eta0 cleaning,
if the longest surviving final-good subsegment is at least as long as the
maximum possible duration of every unprocessed cavity-good candidate, the
search stops. An unprocessed candidate cannot become longer after eta0 cleaning,
so this stopping rule preserves the correct longest-segment choice under the
approved candidate-reset convention.

## FULL CALIBRATION GATE
Full rho_X/rho_Y calibration is allowed only where final_good is true.

## VALIDATION
Synthetic tests must verify:
1. slow cavity and eta0 drifts remain accepted;
2. a reduced cavity block with a deep minimum is rejected even if its average
   remains high;
3. a single deep sample-level eta0 excursion rejects its 0.05 s block;
4. a sustained eta0 collapse longer than the reference window remains rejected
   and cannot redefine the baseline;
5. recovery after a sustained eta0 collapse is accepted;
6. stable eta0 below 0.019 can pass when normal relative to its reference;
7. references remain continuous across raw reader chunk boundaries;
8. the final mask gates rho_X/rho_Y calibration.

The researcher-approved real-data diagnostic over 42 h found approximately
288.40 s of eta0-only rejection with the 80% protected cut. A healthy 3 h file
(0012) remained fully accepted, while the long failure in file 0003 was rejected
for approximately 286 s.

## STATUS
APPROVED BY RESEARCHER AND IMPLEMENTATION AUTHORIZED
