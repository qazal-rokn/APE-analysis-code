# Analysis Contract: Staged APE Cleaning

## REQUESTED
Reject only strong downward cavity-transmission failures and eta0 failures while
preserving slow drift, before the final rho_X/rho_Y calibration.

## CURRENT BEHAVIOR
Raw Dewesoft data are read in bounded synchronized chunks. The canonical
calibration equations are implemented, but production real-data calibration is
gated pending cleaning.

## UNDERSTANDING
Cleaning is intentionally one-sided. Slow drift is not removed. Bad samples are
masked at their original indices and are never deleted and concatenated.

## PLANNED IMPLEMENTATION
1. Summarize cavity transmission in 0.05 s blocks using the block median.
2. Build a centered 60 s local median reference from those block medians.
3. Mark a raw sample cavity-good when
   V_trans >= 0.90 * local cavity reference.
4. Compute I0 and preliminary eta0 only for cavity-good samples.
5. Summarize preliminary eta0 in the same 0.05 s blocks using cavity-good
   samples only.
6. Build a centered 60 s local median eta0 reference.
7. Mark eta0-good when
   eta0 >= 0.90 * local eta0 reference.
8. final_good = cavity_good AND eta0_good.
9. Full rho_X/rho_Y calibration is allowed only where final_good is true.

## SCIENTIFIC DEFINITIONS AND UNITS
The cleaning reference is a local median, not a global mean and not the legacy
90th-percentile reference. Both cuts are lower-side cuts only.

The preliminary eta0 uses the separately approved equation:
eta0 = sqrt(2*sqrt(2)*R2f/(H2*q*G*I0)).

## RISKS AND CONTROLS
A per-chunk local median would create artificial boundary behavior. Therefore
production reference construction is full-record at block resolution while raw
data remain chunked. Only the small block-summary arrays are retained globally.

## WILL NOT DO
No detrending, high-pass filtering, smoothing of science channels, upper-side
cut, interpolation, resampling, gap filling, sample deletion, X/Y combination,
or rho calculation during cleaning.

## FILES THAT WOULD CHANGE
- src/ape_c_lpsd/cleaning.py
- src/ape_c_lpsd/calibration.py
- src/ape_c_lpsd/ingestion/config.py
- configs/dxd_ingestion.example.yaml
- tests/test_cleaning.py
- tests/test_calibration.py
- tests/test_ingestion_config.py
- documentation

## VALIDATION PLAN
Synthetic data contain slow cavity and eta0 drifts plus sharp injected downward
glitches. The slow drift must remain accepted while the injected failures are
rejected. A separate test splits the same synthetic record into chunks to verify
that the full-record reference behaves consistently across chunk boundaries.

## REAL-DATA GATE
Synthetic cleaning tests must pass before this mask is applied to real APE data.

## STATUS
APPROVED BY RESEARCHER AND IMPLEMENTATION AUTHORIZED
