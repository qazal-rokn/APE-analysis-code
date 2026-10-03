# C Method-1 Synthetic Decision Report

Run: `results/synthetic-method1-comparisons/run-f0068be770234bc8a6bf6be6d8f960b9`.
This run used generated inputs only; no APE raw/cache/real-export data was read.

Both methods used the same executable (SHA-256 `fc8a912fba75aee0f26806658ace5414b01ac8887c6fe2675f657b0bf10030ef`), 500 samples at 512 Hz, 24 C nodes in each of 64--128, 128--192, and 192--240 Hz, and matched `-w -2 -p 200 -r 1 -l 50 -a 100 -m 1 -x 1 -M 20` options. Method 0 was direct LPSD. Method 1 used its documented FFT approximation and `-E 10` percent.

All nine tone/noise/mixed cases returned identical frequency nodes. Tone and mixed peaks occurred at their injected C nodes. This is node/peak-location agreement only, not PSD equivalence.

All direct outputs passed strict validation. All method-1 outputs failed it: their final requested row had `AVG=0` and PSD zero. Those rows are preserved in the raw text and run record, not dropped or repaired.

Method-1/direct median PSD relative differences ranged 8.1--51.9 percent. Maximum differences ranged 100 percent to 27,581, largely including low-power/tail bins and the invalid terminal row. White-noise peak locations disagreed in two of three regions.

Direct calls totalled 0.4051 s (0.0450 s mean); method 1 totalled 0.4366 s (0.0485 s mean). Method 1 was 7.8 percent slower in this functional microbenchmark. The child-RSS high-water API reached 11,124,736 KB, but does not isolate per-call memory and cannot compare methods' memory requirements.

An exploratory `-E 1` sparse-grid case did not advance its first approximation block and was terminated; it produced no usable output.

## Decision

Method 1 is not validated or approved for real APE data. It is a distinct FFT-approximation estimator, has a strict terminal-row validity failure, and showed no speed benefit here. A future approval would require a source-level resolution of the terminal row without changing the grid, a non-hanging epsilon/grid domain, larger normalization and noise/tone validation, and a Viper resource benchmark demonstrating practical speedup at equal coverage.
