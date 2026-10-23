# PS-enabled C-LPSD synthetic integration, milestone 1

## Scope
Approved on 23 October 2026: build one first-draft synthetic APE pipeline, with short essential checkpoints; use existing `development/likelihood` and keep `main` immutable until verification.

## Implementation
Existing `spectral.constantq.plan_commands` and original compiled-C PSD mathematics are left unchanged. `spectral.ps_optin.plan_ps_optin` takes an existing six-column plan, an operator-supplied executable, and its independent exact SHA-256. It deep copies the plan, replaces only the C executable, selects `FORMAT "fDSNbri"`, and adds `.raw7` output suffixes to method-2 batches.

`split_ps_batch` validates native C seven columns `f PSD PS Navg Nkernel Re Im` and writes two independent products, preserving original numeric tokens:
- `f PSD Navg Nkernel Re Im` (original six-column PSD contract)
- `f PS` (new direct signal power in rad² if calibrated time series inputs are radians)

No Python FFT/PSD/PS estimator is added. No C execution or cluster submission occurs in the bridge. Each target must be absent; the streaming converter writes exclusive partials and uses hard-links for no-clobber installation. The filesystem must support same-directory hard links; lossless decimal-token preservation is intentional. It does not yet provide the full sixteen-batch join or the run orchestration.

## Essential validation
Three synthetic adapter tests: planner immutability plus digest refusal, valid format conversion with exact token preservation and output overwrite refusal, and non-monotonic frequency rejection with no installed outputs. These unit tests are not a real C execution or science validation.

## Pending for complete integrated synthetic run
1. Retrieve the **validated** opt-in compiled-C source/executable/archive from the archived Mac/VIPER location. Verify archive hash `0297b50652b06be90257e7471dc2cfce84cd15ed45add3129fdd34e0c2d6d48e`, then separately calculate the executable's own SHA-256.
2. Integrate the actual sixteen-batch converter and shared-FFT local/Slurm orchestration, without modifying the frozen C source; run a short X and Y synthetic injection using compiled C.
3. Supply an independently validated physical halo line shape, actual kernel frequency response, and X/Y signal mapping. The provisional QWP peak/mean-square mapping and unit conversion must be tested with injected polarization signals.
4. Call the original vendor background fitting/skew-normal likelihood with consistent PS units and verified normalization; produce likelihood profiles and an explicitly labelled provisional confidence limit. Verify coverage before real-data interpretation.
5. Add an opt-in command-line entry point, parameter manifest, checksums, plots and JSON output. Keep README synchronized at each step.

## Scientific gates
Measured finesse, tau, tau1, rho_DM, C_unit, reference coupling, spectral weights, window responses, quadrature fractions and any axion linewidth assumptions must come from user configuration. **No additional sqrt(f/Q_a)** is inserted in the forward optical response. A fitted PSD [rad²/Hz] must not be mixed with PS [rad²]. Correlated overlapping C-LPSD bins require coverage testing before physical limits.

## Files changed this milestone
`src/ape_c_lpsd/spectral/ps_optin.py`, `tests/test_ps_optin.py`, README and this documentation. Separately, the real-vendor comparison call now uses the identical template in both calls. No changes to vendor or original C source.

## Status
Isolated synthetic adapter tests passed locally. The exact committed GitHub version and the external C runtime have not yet been exercised in this session. The first complete end-to-end run is **not complete**.
