# Synthetic C Grid-Density Feasibility Study

Synthetic input: 128 samples at 128 Hz; Kaiser `-w -2`, PSLL 200 dB, 50%
overlap, `-a 100`, `-m 1`, linear regression, `fmin=9 Hz`, `fmax=36 Hz`.
No APE cache was opened.

| Requested nodes | Result |
|---:|---|
| 4 | finite PSD rows |
| 6 | finite PSD rows |
| 8 | finite PSD rows |
| 10 | finite PSD rows |
| 12 | finite PSD rows |
| 16 | first row: `9 Hz, nan PSD, AVG=0, RBW=0.87074822932 Hz` |

The invalid row is written by the C executable; the strict parser correctly
rejects it. At 16 nodes, the first-node resolution implies a DFT/window length
above the 128 available samples, leaving zero usable estimates. This is a C
estimator feasibility limitation, not a wrapper or parsing error. No invalid row
was discarded, replaced, or accepted.

Smallest tested valid configuration: 4 nodes. Highest tested valid configuration:
12 nodes. The next scientifically honest proposal is a preflight that rejects a
requested C node whose DFT length exceeds the selected uninterrupted segment; it
must report the conflict rather than reduce density. That requires approval.
