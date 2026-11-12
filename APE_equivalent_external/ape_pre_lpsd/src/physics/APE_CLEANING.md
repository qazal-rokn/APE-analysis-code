# Approved APE cleaning definitions

Status: researcher approved.

## Reduced records

The synchronized raw data are partitioned into non-overlapping 0.05 s records
with midpoint times (t_k).

## Cavity stage

For cavity block (k),

[
\bar V_k = \operatorname{mean}(V_{\rm trans}),
\qquad
V_{\min,k} = \min(V_{\rm trans}).
]

The centered 60 s cavity reference is

[
V_{{\rm ref},k}
=
Q_{0.90}
\left(
\{\bar V_j: |t_j-t_k|\le30\,{\rm s}\}
\right).
]

The cavity decision is

[
M_{{\rm cav},k}
=
\left[
V_{\min,k}\ge0.90V_{{\rm ref},k}
\right].
]

The reference is recomputed from provisionally good cavity blocks three times.

## Sample-level eta0

For each raw sample (n),

[
I_0[n]=C_{\rm trans}V_{\rm trans}[n],
]

[
R_{2f}[n]
=
\sqrt{X_{2f}[n]^2+Y_{2f}[n]^2},
]

[
\eta_0[n]
=
\sqrt{
\frac{2\sqrt{2}\,R_{2f}[n]}
{H_2qGI_0[n]}
}.
]

Inside each 0.05 s cavity-good candidate block,

[
\bar\eta_{0,k}
=
\operatorname{mean}(\eta_0[n])_{n\in k},
\qquad
\eta_{0,\min,k}
=
\min(\eta_0[n])_{n\in k}.
]

## Protected eta0 reference

Within one cavity-good segment,

[
\eta_{{\rm ref},k}
=
Q_{0.90}
\left(
\left\{
\bar\eta_{0,j}:
t_k-60\,{\rm s}\le t_j<t_k,\
M_{\eta,j}=1
\right\}
\right).
]

The current block never contributes to its own reference. The eta0 decision is

[
M_{\eta,k}
=
\left[
\eta_{0,\min,k}
\ge
0.80\eta_{{\rm ref},k}
\right].
]

Rejected eta0 blocks do not update the reference. If the accepted-history
window becomes empty during a sustained eta0 failure, the last valid healthy
reference is held until recovery.

At the beginning of every cavity-good segment, eta0 reference state is reset
and the first valid eta0 block mean bootstraps the reference. File boundaries
do not cause a reset when the same cavity-good segment continues across files.

The final mask is

[
M_{{\rm final},k}
=
M_{{\rm cav},k}\land M_{\eta,k}.
]

## Longest-segment search

The full acquisition is first processed with cavity cleaning only.
Cavity-good intervals are ranked by duration. Sample-level eta0 is then
calculated only for candidate intervals that can still exceed the best
surviving final-good interval.

Candidates are processed longest-first. Evaluation can stop when

[
T_{\rm best}
\ge
T_{\rm cavity,next},
]

because eta0 cleaning cannot lengthen an unprocessed cavity-good interval.

No interpolation, smoothing, gap filling, deletion, padding, or concatenation
across rejected intervals is performed.
