# APE C-LPSD / likelihood validation milestone (2026-10-28)

**Status:** frozen local synthetic benchmark, proposed integration into the APE pipeline. This file records the outcome; it does **not** certify a GitHub-hosted end-to-end pipeline or a real-data axion exclusion.

## Scope and non-negotiable scientific rules

- Preserve the upstream `Scalar-Dark-Matter-LPSD` native C estimator and original external vendor likelihood/limit implementation unchanged.
- Follow the 2026 paper's approximation neglecting stochasticity of the dark-matter *signal*: no randomized dark-matter phase or amplitude. Deterministic velocity-broadened expected **power** template.
- Synthetic *noise* residual realizations are permissible for controlled calibration tests.
- `mu=(g/g_ref)^2`, provisional `g_ref=1e-10 GeV^-1`; no additional, unapproved normalization or correction factor.
- Preserve all existing source data, benchmarks, repository history and reference projects.
- Local significances are not look-elsewhere-corrected/global significances.

## Frozen results from local Docker execution

| Component | Outcome / limitation |
| --- | --- |
| Native C-LPSD synthetic spectra | Generated and benchmarked on existing synthetic noise/signal cases; no real APE spectra certified |
| Deterministic axion line-shape template | Original template regenerated with numerical parity; deterministic linewidth model |
| Vendor likelihood signal recovery | 256 residual-noise draws with injected `mu=4`: recovered mean `mu≈4.15`, median `4.0` |
| Local candidate search | Noise-only: 0 selected clusters; deterministic `mu=4`: 7 flagged bins, 1 cluster; peak at `4999.995 Hz`, local `sqrt(q0)=6.7522703365` |
| Vendor `sigma_at_point` versus fixed curvature | 17/17 comparisons successful; maximum relative sigma difference `0.00314991` (0.315%) |
| Vendor median expected 95% limits | 17/17 frequencies. At `5000.000 Hz`: `7.2987030227e-11 GeV^-1`; across 17 frequencies `[4.9869927903e-11, 9.0592789032e-11]` |
| Injected-signal Eq. (19) coverage diagnostic | 256 independent vendor-model noise residual draws, deterministic `mu=4`; 17 exclusions, empirical coverage `93.359375%`, binomial compatibility `p≈0.248`. Not native-C-correlated coverage; none of the draws had negative MLE |
| Positive-MLE observed-limit inversion | All 256 injected-signal trials numerically inverted and 17 exclusions reproduced; native-C synthetic noise at 5 kHz gave `g_UL≈9.0991325945e-11 GeV^-1`. This is **not** a fully validated observed-limit curve |
| Negative-MLE observed-limit convention | **Unresolved.** Signed exact statistic with vendor boundary CDF produced spurious double/no roots; paper signed test and boundary-adjusted convention differ; do not publish 17-frequency observed limits |
| Native-C frequency-bin dependence | Adjacent residual-bin correlations `rho_X≈0.7171`, `rho_Y≈0.7230`; block-resampled versus shuffled score variance ratio `3.6449187235`, std ratio `1.9091670235`. Blocks not independent experiments; these are **not** calibration factors |

## Frozen evidence manifest (local; not uploaded or automatically verified on GitHub)

In the user's Docker workspace `/work`, verified to exist at the time of freeze:

- `ape-paper-eq17-deterministic-linewidth-20261028.json`
- `ape-paper-eq18-local-q0-calibration-20261028.json`
- `ape-paper-eq19-coverage-20261028.json`
- `ape-paper-eq19-observed-ul-20261028.json`
- `ape-paper-vendor-sigma-audit-17freq-20261028.json`
- `ape-validated-local-candidates-expected-ul-20261028.json`
- `ape-native-c-score-correlation-audit-20261028.json`
- `APE_VALIDATION_FREEZE_20261028.md` with SHA-256 digests of the seven evidence files.

The manifest has also been copied to the user's Mac under `~/Desktop/APE_Research/APE_VALIDATION_FREEZE_20261028.md`. This GitHub note is not a substitute for the original data, hashes, source revision or full validation reports.

## Proposed stable local APE entry point

```text
Existing native C-LPSD rho_x/rho_y joined spectra
    -> original vendor background fit + quality status
    -> deterministic APE QWP spectral power template
    -> original vendor likelihood, signed mu estimate, local q0
    -> local candidate clustering (no uncalibrated global p-values)
    -> original vendor sigma_at_point and *median expected* 95% g limit
    -> versioned JSON + CSV, with provenance, warnings and failure status
```

**Integration acceptance:** reuse already approved adapter modules; keep upstream vendor/C source and binaries untouched; run only one end-to-end regression against the frozen 5 kHz benchmark; require exact central expected-limit/candidate parity to recorded precision, manifest consistency, and no production claims. No new statistical experiment is implied by this integration regression.

## Still required for physically meaningful results

1. **Correlated-noise statistical calibration** using representative native C-LPSD observations/independent realizations for candidate global false-alarm control and limit coverage. Resampled contiguous windows are not independent realizations. Define the statistical acceptance target before running.
2. **Physical APE calibration**: confirm complete QWP optical transfer function, cavity response, X/Y quadrature mapping, transmission/PEM calibration, units and coupling normalization, with traceable measured inputs.

Separately decide the **negative signed-MLE observed-limit statistical convention** before enabling an observed-limit curve. Do not silently pick the larger root, alter the original vendor, or substitute a new asymptotic prescription.

## Workflow decision

The frozen synthetic benchmark ends exploratory testing. This proposal records validated components and the fixed integration target. **Do not merge this PR as evidence that the integrated entry point exists until the branch contains and passes the actual integration code and the single frozen-regression check.** Preserve existing `main` contracts for DXD ingestion and any separate scientific approval gates.
