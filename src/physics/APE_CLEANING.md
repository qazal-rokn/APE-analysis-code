# Approved APE cleaning definitions

Status: researcher approved.

The cavity stage follows the earlier validated reduced-block APE method
unchanged. The eta0 stage uses the protected past-only reference validated on
the 42 h APE data set.

## Reduced records

The synchronized raw channels are partitioned into non-overlapping 0.05 s
records with physical midpoint times (t_k).

## Cavity stage

For record (k),

[
\bar V_k = \operatorname{mean}(V_{\rm trans}),
\qquad
V_{\min,k} = \min(V_{\rm trans}).
]

The centered 60 s local cavity reference is

[
V_{{\rm ref},k}
=
Q_{0.90}\left(
\{\bar V_j:\ |t_j-t_k|\le30\,{\rm s}\}
\right).
]

The provisional cavity mask is

[
M_{{\rm cav},k}
=
\left[
V_{\min,k}\ge0.90V_{{\rm ref},k}
\right].
]

The reference is recomputed from provisionally good records three times. The
last refined reference defines the final cavity mask.

## Sample-level eta0

For every raw sample (n),

[
I_0[n]=C_{\rm trans}V_{\rm trans}[n],
]

[
R_{2f}[n]
=
\sqrt{X_{2f}[n]^2+Y_{2f}[n]^2},
]

and

[
\eta_0[n]
=
\sqrt{
\frac{2\sqrt{2}\,R_{2f}[n]}
{H_2qGI_0[n]}
}.
]

For every 0.05 s block (k), define

[
\bar\eta_{0,k}
=
\operatorname{mean}\left(\eta_0[n]\right)_{n\in k},
]

and

[
\eta_{0,\min,k}
=
\min\left(\eta_0[n]\right)_{n\in k}.
]

Only cavity-good blocks can be eta0-good.

## Protected eta0 reference

The eta0 reference is past-only and uses only previously accepted block means:

[
\eta_{{\rm ref},k}
=
Q_{0.90}
\left(
\left\{
\bar\eta_{0,j}:
t_k-60\,{\rm s}\le t_j<t_k, 
M_{{\rm cav},j}=1, 
M_{\eta,j}=1
\right\}
\right).
]

The current block is never included in its own reference.

The eta0 decision is

[
M_{\eta,k}
=
M_{{\rm cav},k}
\land
\left[
\eta_{0,\min,k}
\ge
0.80\,\eta_{{\rm ref},k}
\right].
]

Only accepted blocks update the future eta0 reference. A failed eta0 block is
therefore unable to redefine a long failed state as normal.

If the accepted-history window becomes empty during a sustained failure or
cavity-bad interval, the last valid healthy eta0 reference is held until a
block passes again.

At the beginning of the complete acquisition, the reference is bootstrapped
from the first valid cavity-good eta0 block mean.

## Efficient longest-segment selection

When the analysis goal is specifically to identify the single longest
continuous final-good segment, the full acquisition is first cleaned with the
unchanged cavity criterion. Cavity-good segments are ranked by duration.

Sample-level eta0 cleaning is then applied only to candidate cavity-good
segments that remain capable of exceeding the current best final-good duration.
For this candidate search, the protected eta0 history is reset at the beginning
of each cavity-good candidate and bootstrapped from the first valid eta0 block
mean in that candidate. Within the candidate, the approved past-only q90,
80 percent lower cut, rejected-block exclusion, and held-reference rule are
unchanged.

This reset-at-candidate-start rule is an explicit researcher-approved
optimization for longest-segment selection.

The final mask is

[
M_{{\rm final},k}
=
M_{{\rm cav},k}
\land
M_{\eta,k}.
]

Every raw sample in block (k) inherits these masks. No samples are deleted or
concatenated across rejected blocks. There is no interpolation, smoothing, gap
filling, or upper-side quality cut.
