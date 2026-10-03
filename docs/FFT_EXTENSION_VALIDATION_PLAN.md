# Future Validation Plan: Full-Record FFT Extension

No test in this plan has been run. Synthetic validation and explicit researcher
approval are required before any real voltage sequence is supplied.

| Test | Acceptance criterion |
|---|---|
| One-sided white noise | Across seeded realizations, band-mean PSD agrees with the known variance density within a predeclared Monte-Carlo confidence interval; DC/Nyquist are not doubled. |
| On-bin tone | Peak is at its native FFT bin; integrated one-sided power agrees with injected RMS-squared within a tolerance justified by floating-point precision and window policy. |
| Off-bin tone | Leakage pattern and integrated power agree with the declared window theory; leakage alone is not treated as a failure. |
| PSD/ASD units | `ASD == sqrt(PSD)` at every finite bin and dimensions remain `V^2/Hz` and `V/sqrt(Hz)`. |
| Detrending | `none`, `mean`, and `linear` each produce their documented behavior on deterministic DC/ramp/tone controls. |
| Window normalization and ENBW | Verify `sum(w^2)`, coherent gain, and ENBW against independently calculated values. For `numpy_kaiser_to_be_validated`, compare every sample and normalization value against C `genwin.c` output for multiple N and PSLL values before using any C-compatible wording. |
| C-overlap comparison | With a synthetic uninterrupted voltage input, compare independently saved FFT and C-LPSD products only in their common frequency band. Agreement criteria must be set per estimator/window statistics; equality is not required or assumed. |
| Boundary cases | Test N even/odd, minimum valid N, handoff exactly on/off FFT bin, bins adjacent to DC/Nyquist, and preflight at limit/one-byte-over limit. |
| Channel isolation | Separate X and Y inputs produce separate artifacts, hashes, and provenance; no array combination occurs. |
| Artifact integrity | A second identical requested output refuses overwrite; combined-display metadata references products without creating a merged numerical array. |

Future test reports must state random seeds, NumPy version, configured window and
detrending, sample rate, record length, expected memory, exact acceptance
tolerances, and the C executable/source identity when a C comparison is made.
