# Step 2: Original background fitting on APE X/Y spectra

**Decision: no modification to any vendor source or background algorithm.**

The Step 1 handoff already produces files compatible with original
`scalar_dark_matter_lpsd.LPSDIO.LPSDOutput`:
`APE_<t0>_<t1>_X.h5` and `APE_<t0>_<t1>_Y.h5`.
Pass **each file separately** to the original
`scalar_dark_matter_lpsd.lvk.fit_background` CLI.
The original background fitter, spline model, skew-normal model,
parameter optimization and JSON I/O are reused exactly.

## Run (from APE_equivalent_external/)

```bash
python -m scalar_dark_matter_lpsd.lvk.fit_background \
  --data-path /actual/handoff/APE_T0_T1_X.h5 \
  --output-path /actual/output/ape_background.json \
  --ana-fmin 65 --ana-fmax 10000 --n-processes 1

python -m scalar_dark_matter_lpsd.lvk.fit_background \
  --data-path /actual/handoff/APE_T0_T1_Y.h5 \
  --output-path /actual/output/ape_background.json \
  --ana-fmin 65 --ana-fmax 10000 --n-processes 1
```

Replace T0/T1, paths, and requested frequency interval with measured,
verified values. Both spectra must cover the chosen fitting interval.
Retain original vendor fitting settings unless specifically justified.
Do not pass the entire handoff folder: the original CLI processes
every matching path, including non-spectrum files. Preserve the two
independent records `bkginfo_T0_T1_X` and `bkginfo_T0_T1_Y`.

## Validation performed

In an isolated synthetic experiment, original upstream source code
for `LPSDIO.py`, `lvk/fit_background.py`, `lvk/models.py`,
`lvk/stats.py`, `lvk/hist.py` was used unchanged, with a 750-bin
frequency spectrum for each APE quadrature. Step 1 preserved the
frequency and PSD arrays. The original fitter produced both named
JSON records including `fmin`, `fmax`, `x_knots`, `y_knots`,
`alpha_skew`, `loc_skew`, `sigma_skew`, `chi_sqr`.
Two focused integration tests passed.

**Execution caveat:** the temporary runner lacked `findiff`.
A test-only dummy for its otherwise-unused import was provided to
execute background fitting; no vendor source or scientific function
was patched. The committed regression test does not include that
dummy and requires the normal upstream dependencies.

## Re-run tests on a normal checkout

```bash
python -m unittest discover -s ape_post_lpsd/tests \
  -p 'test_original_background_ape_interface.py' -v
```

No real multi-million-bin APE background fit has been run in this step.
No coupling model, likelihood, X/Y statistical combination,
candidate veto or upper limit has been modified.
