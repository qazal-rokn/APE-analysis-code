# Analysis Contract: Paper-faithful axion-like injection and end-to-end recovery (DRAFT, 2026-10-28)

## REQUESTED
Replace arbitrary seven-bin/correlated-Gaussian demonstrations with one **physical axion-like signal injection** at a prespecified rest frequency, consistent with the 2026 Göttel and Raymond halo / coherence and signal-model assumptions where scientifically transferable to APE. Process real synthetic time-series samples through **unchanged native C-LPSD** and original vendor background, likelihood and inference. Keep it efficient and auditable.

## WHAT WAS WRONG WITH EARLIER TESTS
- The seven-bin Kaiser covariance benchmark tested generic noise dependence, not an axion field and not the paper's line.
- The lab Maxwellian seven-bin likelihood replay tested software but used a coherent-tone pilot and a provisional halo.
- No 2026-paper-faithful, coherence-matched, end-to-end injected axion-like signal has passed.
- Do not treat 12.4% toy false-positive rate as APE false-positive rate.

## PAPER-TO-APE SCIENCE ASSUMPTIONS
- Paper's scalar/LIGO transfer and beta are NOT APE's pseudoscalar QWP response.
- Paper Eq. (5) velocity distribution/virialized halo, escape-speed cutoff, lab/solar motion require specified physical parameters and frame. Do not substitute an arbitrary Gaussian spectral peak.
- Frequency is f = f_rest (1 + v_lab**2/(2 c**2)); rest frequency is LOWER than observed kinetic-broadened frequencies. Characteristic delta f / f ~ 1e-6 is not a hard top-hat full-width.
- The paper's C sampling uses Q=f/delta f (~1e6), native window N_j=fs/(f_{j+1}-f_j), with finite coherent segments and overlap. Actual Q from native log-grid settings; benchmark 5000 Hz => df~0.005 Hz and ~200 s window.
- The paper neglects intrinsic DM signal stochasticity IN ITS LIKELIHOOD on the coherence-matched scale. This is not permission to draw independent power fluctuations in adjacent native C bins; field realization must follow coherent time-domain modes and consistent random phases.
- Original vendor likelihood uses expected signal power plus skew-normal log-power residuals, and its signed optimizer / asymptotic statistics remain unchanged. Statistical adequacy and coverage are outcomes, not guaranteed.
- Real experiment X/Y phase mapping and QWP complex transfer are unvalidated; test must explicitly label any surrogate X/Y phase model.

## PROPOSED SINGLE BOUNDED CAMPAIGN (requires researcher approval)
1. Freeze one physical configuration: f_rest=5000 Hz, Q_C=1e6, lab-halo velocity dispersion, escape cutoff, laboratory motion and annual phase explicitly sourced and serialized; physical mass-frequency conversion and full predicted line distribution; choose deterministic seed. Avoid assumed 7-bin shape.
2. Use a source-derived, truncated, laboratory-boosted Maxwell–Boltzmann velocity field with a *shared* physical realization, consistent phase evolution, and finite coherence. Generate analytic / heterodyned time series plus known broad-band Gaussian noise; write synthetic-only rho_x/rho_y in radians using the existing QWP amplitude envelope, with phase treatment explicitly provisional rather than falsely claiming complete APE physics.
3. Design sampling/recording to satisfy both Nyquist and actual native C window requirements. Avoid claiming that a single 201 s, 0.004975 Hz Fourier grid adequately resolves an entire stochastic halo. Budget length, sample count and memory explicitly. If full resolved-halo synthesis requires >1 coherence window, choose the smallest physically meaningful recording and stop if benchmark budget is insufficient. No enormous blind scan.
4. Execute original PS-enabled lpsd-exec on both channels at exact Q=1e6; native full local window enough for vendor background fit (>=1100 bins), bounded frequency band and memory.
5. Fit original vendor background and evaluate both physically sourced SHM*APE*window expected-power template and original one-bin reference if applicable. No extra K, no arbitrary seven-bin likelihood, no freestanding optical normalization.
6. Predeclare acceptance: frequency localization near actual kinetic-broadened signal support; signal-strength recovery with ensemble/realization uncertainty clearly stated; null control; provenance; measured estimator vs analytic forward PS; original vendor numeric parity. A single halo realization is NOT a coverage proof; do not claim global detection significance or physical 95% limit.
7. One end-to-end run with a small number of predeclared stochastic field realizations and null controls, reusing native C FFTs wherever valid. No repeated source normalization tests. Preflight before runtime; cap time/memory.
8. Write one report and archive run-specific outputs into a fresh no-overwrite path.

## FILES AFFECTED IF APPROVED
A new isolated physically sourced halo-time-series injection/preflight tool under tools/ or src/physics/, a focused integration acceptance script and test(s), one config, and synchronized README / docs/EQUATION_REFERENCE.tex / this contract. Existing C, vendor, calibration and successful artifacts remain immutable.

## IMPORTANT UNRESOLVED PHYSICS
Whether and how the 2026 scalar paper's phase/coherence construction transfers to pseudoscalar APE, the actual optical QWP X/Y transfer phase, exact halo parameter sources, live vs frozen geometry, noise amplitude and justified injection strength, and what constitutes a valid frequency- and coupling-recovery tolerance. These must be explicitly approved or treated as declared benchmark choices, never as measured APE facts.

## REAL DATA
Not authorized.

## STATUS
DRAFT. Researcher requested correct experiment; scientific implementation not started pending approval of assumptions and contract. No C or other pre-existing file modified.
