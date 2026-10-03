# Units and normalization

The bridge transports cached float64 `rho_x` or `rho_y` samples unchanged with
`-x 1`. It does not recalibrate or renormalize. C output `D` is preserved as the
one-sided PSD; ASD is only `sqrt(PSD)`. If the input unit is `U`, units are
`U^2/Hz` and `U/sqrt(Hz)`.

The approved C backend is intentionally not numerically identical to the prior
Python one-prefix estimator: it uses C LPSD averaging and `-r 1` linear
regression, whereas the previous Python adaptation used mean-only detrending.
Calibration, mask, selected segment, units, and requested logarithmic density
are preserved; the spectral estimator is the intended change.

`1/T=fs/N` remains the diagnostic one-cycle frequency. It is not necessarily a
valid C-LPSD node. The wrapper obtains the actual C-reported Kaiser `SBIN` using
a synthetic C probe with the configured window, PSLL, overlap, and averaging,
then uses `fmin = SBIN/T * [1+(1-overlap/100)(m-1)]`. With `m=1`, this is
`SBIN/T`. This is a C-estimator feasibility limit from finite window width, not a
cleaning cut or missing-data rule. `-a 100` remains desired; realized averages
are output in column `N` and must be recorded.
