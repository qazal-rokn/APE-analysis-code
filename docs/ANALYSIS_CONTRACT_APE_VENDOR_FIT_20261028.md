# APE single-bin vendor MLE to expected UL (provisional; 28 Oct 2026)

## REQUESTED
Connect the original program's likelihood optimizer, one-sided q0 and 95% expected
UL in APE development branch without changing original vendor statistics.

## CURRENT BEHAVIOR
The original likelihood and expected-UL adapter passed direct tests in Docker;
the vendor max_lkl_executable.calc_max_lkl passed a one-bin synthetic test.
The APE controller previously used a finite positive grid and did not report
vendor fitted sigma_mu or expected limits.

## UNDERSTANDING
Use original signed-mu maximization, Fisher-derived fitted sigma_mu,
original q0=2(Lmax-L0) if muhat>0 else zero, and original expected-UL function.
APE parameter is mu=(g/g_ref)^2. Preserve its g_ref scaling exactly.

## PLANNED IMPLEMENTATION
Add one new APE adapter that calls calc_max_lkl unchanged, passes sigma_mu to
existing unchanged-vendor expected-UL adapter, returns q0 and provenance.
Restrict pilot to the original default single-bin peak_shape=[1]. Require
explicit optimizer scan settings. Do NOT alter existing candidate scanner.

## SCIENTIFIC DEFINITIONS AND UNITS
Y, bkg are log power; peak_norm*peak_shape is signal PS rad^2 at g_ref;
mu is dimensionless, sigma_mu is the optimizer uncertainty; g_UL is GeV^-1.
The expected 95% limit is NOT an observed limit and does not measure 95%
coverage for APE physical stochastic axion signals.

## RISKS AND CONTROLS
Vendor get_valid_mu_ranges does not multiply peak_shape in its range scan,
so multi-bin use requires an explicit scientific validation, not silent support.
Vendor optimizer may fail/return NaNs at zero or near a boundary; fail closed.
Correlated C bins, original signal shape provenance and physical calibration
remain unvalidated here. The test is single-bin software parity only.

## WILL NOT DO
No vendor/C edits; no raw data; no native C run; no background refit;
no modifications to original APE calibration or scientific source definitions;
no Git branch switching or automatic commit/push.

## FILES THAT WOULD CHANGE
New ape_vendor_fit.py, test_ape_vendor_fit.py, this analysis contract,
append a short scope note in README and add q0/UL mapping in EQUATION_REFERENCE.tex.

## VALIDATION PLAN
Numerical parity to original vendor MLE on earlier deterministic fixture,
original vendor expected-limit parity, and refuse multibin inputs.
Run only existing software tests. Scientific coverage remains separate.

## UNRESOLVED DECISIONS
Physical APE multi-bin halo/PS-kernel template validation, candidate global
trial factor calibration, noise correlations, and observed UL definition.

## REAL-DATA GATE
Not authorized.

## STATUS
Provisional integration under existing researcher direction; not merged or
approved for physical inference.
