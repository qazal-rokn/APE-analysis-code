# APE vendor MLE to local candidate report (single-bin pilot, 28 Oct 2026)

## REQUESTED
Use the original paper's unchanged maximum-likelihood fitting, one-sided q0,
expected upper limits and APE candidate search in the development branch.

## CURRENT BEHAVIOR
A true vendor MLE and vendor expected-UL adapter passed numerical parity in
Docker. Existing APE candidate preselection implements the vendor q0 and fixes
the original vendor's missing terminal-bin cluster. Candidate/UL bridge absent.

## UNDERSTANDING
The existing single-bin vendor optimizer is the only approved pilot.
Use the original vendor likelihood/fit and vendor expected-limit method.
No invented signal statistic, global trial factor, cross-quadrature veto,
physical halo template or external one-bin response correction.

## PLANNED IMPLEMENTATION
One thin report function taking multiple already-prepared single-bin candidate
windows and their exact native-grid center frequencies; fit each using the
previously tested vendor adapter and call existing select_candidates. No
changes to either vendor or candidate statistic. Return provisional expected
upper limits alongside local candidate q0 values. No file I/O or full search.

## SCIENTIFIC DEFINITIONS AND UNITS
mu=(g/g_ref)^2, q0=2*(Lmax-L0) for muhat>0 and q0=0 otherwise.
The user supplies a positive local sigma threshold. sqrt(q0) is a LOCAL
ranking statistic only, not a physical detection significance. The original
vendor upper limit is asymptotic EXPECTED, not observed; g_UL is GeV^-1.

## RISKS AND CONTROLS
Correlated C bins, full trial distribution, halo PS shape and QWP X/Y response
are unvalidated here. A single-bin default is not a realistic multi-bin
axion injection. No 95% APE physical coverage claim.

## WILL NOT DO
No real APE data, native C, vendor edits, branch changes, previous-result
overwrite, upper-limit method replacement or arbitrary signal normalizations.
No changes to existing finite-grid controller.

## FILES THAT WOULD CHANGE
New report adapter, one focused test, this contract, README and equation notes.

## VALIDATION PLAN
Local ranking/terminal cluster unit check, direct vendor-fit parity on the
previous deterministic fixture, invalid-input refusal, clean Git review.

## UNRESOLVED DECISIONS
Multi-bin original-paper signal shape versus native C power kernel, actual
background fitting for selected bins, effective trial significance, physical
injection recovery and frequentist coverage.

## REAL-DATA GATE
Not authorized.

## STATUS
Provisional single-bin integration checkpoint, not merged into main.
