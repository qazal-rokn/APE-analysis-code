# APE vendor expected upper-limit adapter (28 October 2026)

## REQUESTED
Reuse the original external paper's statistical functions unchanged and test
candidate searching and 95% upper limits before integration with full APE.

## CURRENT BEHAVIOR
APE already delegates vendor log-likelihood and has a candidate q0 preselector.
Original vendor get_upper_limits_approx has been tested directly in Docker.
No main APE controller integration or observed UL calculation is authorized here.

## UNDERSTANDING
Preserve the vendor asymptotic expected limit, distinguish it from observed
profile-likelihood inversion and from physical 95% coverage.

## PLANNED IMPLEMENTATION
Add exactly one guard-checked adapter that calls vendor
get_upper_limits_approx(sigma_mu, alpha_CL, tol), and tests its return
convention using the previously verified numerical result. No new optimizer.

## SCIENTIFIC DEFINITIONS AND UNITS
APE defines mu=(g/g_ref)^2, sigma_mu is the UNRESCALED fitted uncertainty
of mu in that exact convention. The vendor returns sqrt(mu_UL) and its
uncertainty; return g_UL=g_ref*sqrt(mu_UL) in GeV^-1. The adapter DOES NOT
estimate sigma_mu: the caller must supply independently validated fit output.

## RISKS AND CONTROLS
The original vendor result is an asymptotic expected limit, not observed UL.
No 95% coverage for physically simulated APE axions is implied. A finite-grid
candidate optimizer does not replace vendor fitted mu/sigma; existing APE
candidate implementation and original vendor clustering remain unchanged.

## WILL NOT DO
No edits to original vendor/C sources, calibration, signal template, native PS,
main controller, candidate selection, raw data, or earlier results.

## FILES THAT WOULD CHANGE
New adapter module, new focused test and this contract only.

## VALIDATION PLAN
Test direct parity with the original vendor's known sigma_mu=1 result,
parameter scaling, invalid input refusal, Git/source immutability.

## UNRESOLVED DECISIONS
How original vendor maximum-likelihood sigma_mu will be obtained for the
physical APE template, and how native-C correlations affect coverage.

## REAL-DATA GATE
Real-data upper limits not authorized.

## STATUS
Narrow integration approved in ongoing researcher instruction; provisional.
