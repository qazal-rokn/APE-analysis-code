# C-LPSD Metadata-Only Feasibility Map

Inputs: selected segment length $N=67,199,997$, $f_s=20,000$ Hz, C-probed
Kaiser SBIN=8.15, PSLL 200 dB, 50% overlap, and minimum averages 1. No signal
array was opened.

The C node planner uses $L=\mathrm{round}(f_s/\delta f)$. Thus a node needs
$L\le N$, or locally $\delta f/f\ge f_s/(Nf)$. The corresponding independent
local maximum relative density is approximately $Q_{max}=Nf/f_s$.

| Frequency (Hz) | local $Q_{max}$ |
|---:|---:|
| C lower bound 0.0024256 | 8.15 |
| 0.01 | 33.6 |
| 0.1 | 336 |
| 1 | 3,360 |
| 10 | 33,600 |
| 100 | 336,000 |
| 1,000 | 3,360,000 |

For an individually chosen local relative resolution $Q=10^6$, feasibility
starts at $f=f_s Q/N=0.297619$ Hz. For the currently requested whole-range C planner with $J_{des}=17,330,046$,
there is no lower-frequency choice that makes every node from that lower limit to
Nyquist feasible: the smallest required window remains above $N$ (its best case
is roughly $9.4\times10^7$ samples). Therefore there is no valid answer to
"lowest frequency at which the current whole-range Q=1e6 grid is feasible."

The C program can naturally cover the full C-valid range only with a much
coarser node count (about 122 nodes by the strict finite-segment bound at the
lowest C frequency). This is an assessment, not an automatic recommendation.

Beginner interpretation: fine frequency spacing needs long chunks of data. At
very low frequency this C implementation asks for more samples than exist in the
longest uninterrupted record. The data are not missing; the requested estimator
grid is too fine for its finite-window rule. A researcher must choose a new,
explicitly approved scientific path before a spectrum is attempted.
