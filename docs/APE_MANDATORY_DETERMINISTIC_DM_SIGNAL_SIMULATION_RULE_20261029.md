# Mandatory 2026-paper synthetic-signal assumption

This is a standing rule for APE paper-based likelihood/injection validation.

1. The dark-matter **signal** is represented by a deterministic, velocity-broadened **expected power spectrum**, using the established truncated-Gaussian halo distribution, optical QWP response, and appropriate C-LPSD window-power response.
2. Do **not** randomize axion phase or amplitude, do not generate independent random axion coherence-window realizations, and do not present a single coherent sinusoid as the paper-based physical axion injection.
3. Random detector-noise realizations are permitted, and matched noise-only and injected tests should share the same noise realization.
4. The power-domain signal model alone does not uniquely determine a time series. Generating a field-level time series would need an **additional explicit phase convention** inconsistent with the stipulated no-random-phase requirement if phases are drawn randomly. Keep a separate C-LPSD software smoke test rather than claiming a power injection exercises C.
5. The original compiled C estimator, vendor background/likelihood/upper-limit code, accepted cleaning and calibration, and prior frozen data remain unchanged. Never insert an ad hoc `f/Q_a` or `N^-1/4` factor into the native PS-to-likelihood path.
6. The power-spectrum fixture is a scientific **diagnostic**, not a validated physical APE detection or exclusion. Any updated optical or statistical convention needs an explicit separate scientific review.

Historical note: the prior 250-Hz coherent-tone package was only a software smoke fixture and is **not approved** for paper-based axion-signal validation. This project rule supersedes its interpretation, without modifying or deleting any earlier benchmark files.

## Mandatory pre-simulation paper review and statistical scope

Before **each new APE axion-related simulation**, re-check the cited 2026 paper (DOI 10.1038/s41598-025-33428-2), its relevant equations and approximations, and the original vendor code. Record which assumptions are source-defined and which numerical choices are synthetic. Recheck that the dark-matter expected **power** is deterministic and that no random dark-matter phase, amplitude, or coherence-window draws are introduced. Random detector noise is permitted.

Do not conflate the original vendor's **default one-frequency-bin profile likelihood** (`peak_shape=np.ones(1)`) with an APE **multi-bin halo-template** likelihood. Correlations *between tested frequency hypotheses* do not alone invalidate a correctly specified local one-bin likelihood; correlations *between observations included in one joint likelihood* may. Independent local p-values, global candidate trials, and physical observed upper-limit coverage require distinct calibration. Do not claim full APE or native-C validation from a power-level statistical test. The final single-bin Monte Carlo is a declared reconstruction, not an exact reproduction of unreported paper Monte Carlo parameters.
