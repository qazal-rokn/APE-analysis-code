# Approved APE cleaning definitions

Status: researcher approved.

The cleaning architecture follows the earlier validated APE reduced-block
cavity-quality method, with a researcher-approved relative eta0 rule replacing
the historical fixed eta0 >= 0.019 threshold.

## Reduced records

The synchronized raw channels are partitioned into non-overlapping 0.05 s
records. Each record carries its physical midpoint time.

## Cavity stage

For record k define

\[
\bar V_k = \operatorname{mean}(V_{\rm trans})
\]

and

\[
V_{\min,k} = \min(V_{\rm trans}).
\]

The local cavity reference is the centered 60 s 90th percentile of eligible
record averages:

\[
V_{\rm ref,k}
=
Q_{0.90}\left(
\{\bar V_j:\ |t_j-t_k|\le 30\,{\rm s}\}
\right).
\]

The provisional cavity mask is

\[
M_{{\rm cav},k}
=
\left[
V_{\min,k}\ge 0.90 V_{\rm ref,k}
\right].
\]

The reference is recomputed from provisionally good records three times. The
last refined reference defines the final cavity mask.

## Preliminary eta0

For each reduced record,

\[
\bar R_{2f,k}
=
\sqrt{\bar X_{2f,k}^{,2}+\bar Y_{2f,k}^{,2}},
\]

\[
I_{0,k}=C_{\rm trans}\bar V_k,
\]

and

\[
\eta_{0,k}
=
\sqrt{
\frac{2\sqrt{2}\,\bar R_{2f,k}}
{H_2 q G I_{0,k}}
}.
\]

Only cavity-good records are eligible to define the eta0 reference.

## Eta0 stage

The eta0 reference is

\[
\eta_{{\rm ref},k}
=
Q_{0.90}\left(
\{\eta_{0,j}:\ |t_j-t_k|\le30\,{\rm s},\ M_{{\rm cav},j}=1\}
\right).
\]

The provisional eta0 mask is

\[
M_{\eta,k}
=
M_{{\rm cav},k}
\land
\left[
\eta_{0,k}\ge0.90\eta_{{\rm ref},k}
\right].
\]

The eta0 reference is recomputed from provisionally eta0-good and cavity-good
records three times. The final mask is

\[
M_{{\rm final},k}=M_{{\rm cav},k}\land M_{\eta,k}.
\]

Every raw sample in record k inherits these masks. No samples are deleted or
concatenated across rejected records. There is no interpolation, smoothing,
gap filling, fallback reference, or upper-side quality cut.
