# APE signal and likelihood connection, draft 2026-10-23

## Approved scope
Connect the verified C PS executable interface, an explicitly parameterized axion spectrum, and original vendor likelihood on development/likelihood only. Defer synthetic-data generation and end-to-end execution. Update README in the same change. No changes to main, original vendor, compiled C, data, or historical outputs.

## Physics
The existing nominal QWP model is
R(f) = tau * C_unit * sqrt(2*rho_DM / ((pi/F_s)^2 + 4 sin^2(pi f tau1))).
For the existing provisional peak-amplitude convention:
T_qj = (g_ref^2/2)*alpha_q*sum_k h_k*R(f_k)^2*W_qjk [rad^2].
The vendor model uses B_qj + mu N_qj S_j with N_qj S_j = T_qj and mu=(g/g_ref)^2. Do not insert sqrt(f/Q_a). All QWP and axion parameters are explicitly user supplied.

## Halo
The optional unboosted Maxwellian benchmark implements normalized speed density p(v) proportional to v^2 exp(-v^2/v0^2) with v^2 = 2*c^2*(f/f_rest-1). It computes discrete bin probability masses using the analytic CDF and refuses insufficient frequency coverage. This excludes lab velocity, escape truncation and modulation. A physically approved laboratory-frame halo must replace it for actual limits. Supplying a separately validated halo probability array remains supported.

## C boundary
The PS opt-in adapter verifies the exact supplied C executable hash. The archived tarball digest is NOT the binary digest. In-session execution of this binary cannot be claimed. The actual C-LPSD power kernel W_qjk is not yet derivable from the frozen backend and must be supplied/validated; never replace it with a NumPy FFT/Welch estimate.

## Statistics and safety
The APE glue makes arrays (M,N) and calls the unchanged vendor likelihood; no vendor statistical algorithms were reimplemented. Background PS and observed PS must both be rad^2 and strictly positive. X and Y fractions, correlated bins, peak/RMS normalization, finite-coherence line shape and inferential coverage remain unverified. No physical confidence bounds will be reported until calibration and injection tests.

## Tests and status
Only minimal static contract tests have been added for the halo coverage and APE likelihood normalization mapping. No genuine C or vendor runtime and no synthetic end-to-end data generation performed in this milestone. README updated. Main unchanged.
