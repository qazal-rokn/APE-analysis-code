# Mandatory APE simulation contract: 2026 paper signal approximation

For any APE synthetic injection, recovery, end-to-end statistical test, and reported sensitivity comparison:

1. **Neglect stochasticity of the dark-matter signal**, as in the approved 2026-paper analysis approximation. No random axion phases, no random axion amplitude realizations, no independent random phase per coherence time, no random discrete halo-mode sampling used as the signal template.
2. Use the existing **deterministic, velocity-broadened expected axion SIGNAL POWER** template: fixed truncated isotropic Galactic 3D Gaussian benchmark and lab boost, convolved with the existing approved native C-LPSD Kaiser power-window kernel. Use only the existing provisional QWP amplitude response, with its provisional status explicit.
3. **Random detector noise is allowed**, including repeated independent noise-only trials for null statistics and coverage assessment. It is the NOISE realization, not the DM signal power, that varies between trials.
4. Power-domain injection must be explicit: `P_observed(f) = P_noise(f) + (g_injected/g_reference)^2 * T_ape(f; g_reference)`, with native PS [rad²], where `T_ape` already incorporates the existing optical response, halo shape, window, and X/Y template. No additional `f/Q_a`, `sqrt(f/Q_a)`, or ad hoc `N^(-1/4)`.
5. Keep real raw data, cleaning, calibration, original C, original external vendor likelihood/limits, and frozen benchmarks **unchanged**. New test fixtures must have separate paths and explicit hashes.
6. A **deterministic power-domain injection after C-LPSD** tests the vendor statistical chain without changing the C noise estimator. It is **not** a time-domain axion realization, and by construction does not validate physical propagation of an axion signal through C. Never claim an end-to-end injected time-series test from such a fixture.
7. A coherent sinusoid or the older `finite_random_phase_halo_plane_wave_superposition` is **not an approved paper-consistent signal fixture** for the main validation. Such tests must be labeled as distinct diagnostics, not substituted.
8. The existing X/Y response and QWP normalization remain provisional until separately validated, and the global trials, correlated C spectral-bin coverage, and observed-limit conventions remain open.
9. Every synthetic signal manifest must explicitly include `random_DM_phase=false`, `random_DM_amplitude=false`, `stochasticity_of_DM_signal_neglected=true`, `injection_domain`, `source_paper_doi`, and source code/configuration hashes.

This contract is persistent in the repository so a future test cannot silently revert to random dark-matter phases. It does not modify or supersede scientific code.
