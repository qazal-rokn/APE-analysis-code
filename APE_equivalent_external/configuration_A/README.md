# Configuration A — minimal original-vendor analysis

**One working branch:** `experiment/ape-configuration-a-summed-quadratures-20261119`.

Configuration A has only two Python modules:

- `combine_quadrature_psd.py`: the necessary APE-only handoff that computes
  `PSD_Q = PSD_X + PSD_Y` in **linear** PSD units, then
  `logPSD_Q = log(PSD_Q)`. It preserves the native C grid and produces an
  HDF5 file compatible with the **original** `LPSDOutput`.
- `run.py`: a small CLI that calls the **original**
  `scalar_dark_matter_lpsd.lvk.fit_background` for Q, then passes Q into
  the existing APE-only detector-response adapter
  `ape_post_lpsd.ape_max_likelihood`, which invokes the **original**
  `max_lkl_executable.calc_max_lkl` and `stats.log_likelihood`.

The original vendor `LPSD/` C source, background fitter, likelihood,
optimizer, and upper-limit/statistical modules are unchanged. The model
uses the original `peak_shape=[1]` single-bin convention. No separate
Monte Carlo/halo/window fit is used in Configuration A.

## 1. Combine X/Y and calculate a NEW background using original vendor code

From `APE_equivalent_external/`:

```bash
python -m configuration_A.run background \
  --x-h5 /work/APE_0_15435_X.h5 \
  --y-h5 /work/APE_0_15435_Y.h5 \
  --output-dir /work/combined_Q \
  --ana-fmin 65 --ana-fmax 8000
```

This produces `/work/combined_Q/APE_0_15435_Q.h5` and
`/work/combined_Q/background_Q.json`. Previous X/Y backgrounds are
not valid for Q and are not modified. The vendor fitter uses its
unchanged defaults unless optional `--n-processes` or
`--segment-size` is explicitly supplied. Refit Q before likelihood.

## 2. Run original vendor single-bin likelihood with APE response

Supply a precomputed APE `response.npz` containing
`frequency_hz`, `peak_norm_x`, `peak_norm_y`, and `metadata_json`.
It is **not generated or modified** in Configuration A.
The existing minimal APE response adapter checks the input and sets
`peak_norm_Q = peak_norm_x + peak_norm_y` before calling the original
optimizer. The optical response must be in the same calibrated angle PSD
units as the summed spectrum.

```bash
python -m configuration_A.run likelihood \
  --spectrum /work/combined_Q/APE_0_15435_Q.h5 \
  --background-json /work/combined_Q/background_Q.json \
  --response-npz /work/ape_response.npz \
  --output-npy /work/combined_Q/likelihood_Q.npy \
  --ana-fmin 65 --ana-fmax 8000 \
  --min-log10mu MIN --max-log10mu MAX
```

Choose `MIN` and `MAX` for the response's actual units.
No fabricated response or arbitrary optical constants are supplied.

**State:** code updated on GitHub; no real-data Q background refit or
full-band likelihood execution has been performed here. An actual
calibrated APE response NPZ is required for meaningful coupling results.
Configuration B stays separate; existing exploratory scripts elsewhere
in the repository are not imported by these two modules.
