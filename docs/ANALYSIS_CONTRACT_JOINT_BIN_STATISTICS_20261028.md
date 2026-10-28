# Analysis Contract: APE joint-bin statistical validation (DRAFT, 2026-10-28)

## REQUESTED
Evaluate original paper likelihood assumptions, quantify C-LPSD neighboring-bin and X/Y dependence, and perform ONE consolidated, efficient synthetic inference validation for strength recovery and false-positive behavior. No repeated C normalization tests.

## CURRENT BEHAVIOR
Original vendor stats.py computes residuals r[q,j] = log P[q,j] - log(B[q,j]+mu*T[q,j]), evaluates skew-normal PDFs, and sums individual log densities over q and j. No off-diagonal bin or cross-X/Y covariance is modeled. Original vendor's paper states signal stochasticity is neglected under its coherence-time-matching approximation; this is not an APE-specific physical validation.
APE's physical_multibin_template mode evaluates unchanged original vendor log-likelihood on a prescribed coupling grid for synthetic inputs, without calibrated significance or limits. A 5 kHz coherent tone and the 7x7 Kaiser power response have been checked independently. The historical pilot uses Q_eff ~500149, not Q=10^6.

## INITIAL ANALYTIC RESULT (SOURCE + DERIVATION)
For zero-mean stationary *circular complex Gaussian white noise* and frequency estimators using the same centered window w(t), Cov(|Z_j|^2,|Z_k|^2)/sqrt(Var...)=|<w_j,w_k>|^2 / (<w_j,w_j><w_k,w_k>). This is a theoretical **power** correlation, not covariance of vendor skew-normal **log residuals** and not a measurement of real APE noise.
Using the original C Kaiser window with PSLL=238.13 dB, alpha=9.549175266675151, and 200 s same-length windows on 0.005-Hz grid gives approximate adjacent-lag correlations 0.849375, 0.519776, 0.228234 at 1, 2, 3 bin separations. C-Q variability, window-centering/segment offsets, unequal lengths, detrending, actual noise spectrum and X/Y cross-coupling can change these values. Do not substitute this ideal covariance for measured log-power covariance or into original vendor source.

## UNDERSTANDING
Scientific goal is not to prove the unchanged vendor code executes again. It is to determine how much dependence modifies false alarms, effective information and nominal upper-limit coverage. The physical stochastic axion time-domain and X/Y phase model must be explicitly defined and validated before physical g recovery claims.

## PLANNED IMPLEMENTATION (REQUIRES RESEARCHER APPROVAL)
1. Add one opt-in, **read-only** statistics-validation program outside the original vendor/C and outside existing prior result paths. Reuse the existing precomputed C kernel, original log-likelihood function and source-verified halo/QWP convention.
2. Record the actual baseline dataset conditions and signal-generation prescription, seed and noise representation. If appropriate validated *noise-only* independent realizations are unavailable, use a clearly labeled analytic/stochastic covariance benchmark first; DO NOT claim to have measured APE noise covariance from a single coherent pilot or from seven spectral powers.
3. For a prespecified local candidate window, build one **controlled ensemble** with null and injected signals using shared (rather than independent per-bin) underlying random degrees of freedom. Preserve original expected field coherence assumptions; synthetic benchmark noise model declared explicitly. No new C execution if a statistically legitimate frequency-domain simulation is sufficient for this diagnostic. Do not invent per-bin independent DM powers.
4. Compute: covariance of *vendor log-power residuals* empirically; X/Y covariance or a documented zero-correlation premise; fitted mu/g bias and uncertainty; null q0 distribution and empirical false-positive rates; observed 95% coverage where a defined interval exists, comparing original vendor product likelihood with a separate *diagnostic* dependence-aware comparator. Do not claim global trials calibration from a local test.
5. Predeclare ensemble count, effect size/precision targets, seed, runtime and stop criteria in config. Vectorize and reuse cached quantities; run **one consolidated** acceptance test only. Reject silently truncated/renormalized templates and nondeterministic copies of existing data.
6. Summarize a scientifically binary verdict: original vendor approximation demonstrably adequate in the specified benchmark or not established. Do not silently replace original statistics; any adoption of a correlated likelihood is a separately authorized scientific methodology decision.

## SCIENTIFIC DEFINITIONS AND UNITS
Native C power is rad^2; existing forward T[q,j] = g_ref^2/2 * alpha_q sum_k h_k R_QWP^2(f_k) W[q,j,k]. No extra empirical K normalization. mu=(g/g_ref)^2, signed vendor MLE for internal modeling; physical coupling g>=0. Local q0=0 if muhat<=0, otherwise 2(logLmax-logL0). Empirical null q0 CDF/false-positive control must not be replaced by asymptotic chi-square assumptions without evidence.

## RISKS AND CONTROLS
An independent-bin surrogate ensemble cannot validate actual overlapping C bins. A stationary Gaussian-noise covariance of power alone does not determine skew-normal log residual covariance. LIGO scalar response and independent detectors are not APE QWP and same-detector X/Y; dependence must be separately modeled or measured. Physical SHM input and QWP phase remain provisional. Initial coherent signal-only tests are not null ensembles.

## WILL NOT DO
No new large C jobs, no changes to original C/vendor/source, physics/calibration, frozen caches, main or real data; no observed g limits/detections; no unapproved covariance-model substitution.

## FILES THAT WOULD CHANGE
One isolated opt-in validation tool under tools/, one focused contract test under tests/, one report under docs/ and a short method note. Exact paths and test profile reviewed before implementation.

## VALIDATION PLAN
Single predeclared deterministic small ensemble with comparison to closed-form Gaussian limiting case and vendor-likelihood parity. Run from Linux venv with one command. No repeated integration pilot.

## UNRESOLVED DECISIONS
The physical APE axion stochastic field and X/Y phase covariance prescription, whether actual independent noise-only C spectra exist and can be reused, which coverage criterion to target, and what diagnostic reference correlated likelihood to use.

## REAL-DATA GATE
NOT AUTHORIZED.

## STATUS
DRAFT for researcher approval. Read-only source audit and idealized analytic benchmark completed. No new statistical analysis implementation and no inference validation performed.
