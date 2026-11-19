# Full 15,435-second APE native-grid observed likelihood scan

**Purpose:** use the already completed X/Y four-hour spectra and *original
full-grid vendor background fits* to evaluate the observed axion likelihood
at **every eligible native Constant-Q frequency bin** between 65.1 and
8,000 Hz. This supersedes the 64-point engineering survey for frequency
coverage. It does not rerun the C estimator or background fitting.

**Interpretation:** a provisional *observed composite-likelihood scan*,
NOT a calibrated exclusion limit, significance, candidate claim, or a
publication-ready axion response. The original vendor
`scalar_dark_matter_lpsd.lvk.stats.log_likelihood` is evaluated directly
on the original logPSD, original background spline and original fitted
skew-normal residual parameters. The dense optimizer is
`scipy.optimize.minimize_scalar` in one scaled signed amplitude variable
to avoid the original vendor's expensive seeds and frequent Fisher
`Spacing must be > 0` failures. This changes the optimizer only,
NOT the likelihood function. Comparison against vendor optimizer at selected
points remains an acceptance check.

At every rest-frequency hypothesis the code builds the previously
validated 81-bin analytic Kaiser-PSD kernel, applies the fixed laboratory
halo benchmark and QWP transfer from `physics/qwp_signal.py`, and scans
**separate X/Y conditional unit-power signal projections**. Do not
jointly combine X/Y, or interpret `mu_hat`, q0 as calibrated. In
particular, neighboring native C-LPSD bins overlap and are correlated;
the vendor product PDF assumes independence. The halo Maxwellian
parameterization, experimental optical transfer (including phase and X/Y),
and production PSLL=238.13 dB are not certified.

No bulk `K_D` matrix is stored: each axion mass uses only an 81x64
temporary kernel (by default). A single dimensionless sampled halo is
reused and scaled to each rest frequency, preventing MC jitter between
neighboring hypotheses. The Kaiser normalization quadrature is cached in
each worker. Original HDF5/native-grid verification and background
segmentation are performed for each task chunk, not separately for each mass.
Each task writes one atomic compressed NPZ checkpoint containing results,
validity flags, native indices and provenance metadata. No source data
is changed.

## Run on VIPER

```bash
cd /ptmp/qarokn/APE-analysis-code
git pull --ff-only
cd APE_equivalent_external
PY=/ptmp/qarokn/ape-python-venv-20261023/bin/python
ROOT=/ptmp/qarokn/APE_4p28h_PSD_OriginalVendor_PostLPSD_20261117
export PYTHONPATH="$PWD"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1

"$PY" -m unittest discover -s ape_post_lpsd/tests \
  -p 'test_ape_full_native_likelihood_scan.py' -v

"$PY" -m ape_post_lpsd.ape_full_native_plan \
  --spectrum "$ROOT/handoff/APE_0_15435_X.h5" \
  --task-bins 10000 --fmin 65.1 --fmax 8000
```

First run a **200-mass compute-node benchmark**, separate from
production output, and inspect the elapsed seconds and fit quality:

```bash
"$PY" -m ape_post_lpsd.ape_full_native_likelihood_scan \
  --root "$ROOT" \
  --outdir "$ROOT/diagnostics/full_native_pilot_200" \
  --task-id 0 --task-bins 200 --fmin 65.1 --fmax 8000
```

If that pilot finishes and matches the vendor results/expected numerical
behavior, submit production jobs from a VIPER Slurm submit shell:

```bash
OUT="$ROOT/diagnostics/ape_full_native_4p28h_v1"
mkdir -p "$OUT/logs"
TASKS=$("$PY" -m ape_post_lpsd.ape_full_native_plan \
  --spectrum "$ROOT/handoff/APE_0_15435_X.h5" \
  --task-bins 10000 --fmin 65.1 --fmax 8000 |
  "$PY" -c 'import json,sys;print(json.load(sys.stdin)["tasks"])')

sbatch --array="0-$((TASKS-1))%8" \
  --export="ALL,ROOT=$ROOT,OUT=$OUT" \
  --output="$OUT/logs/%A_%a.out" \
  ape_post_lpsd/ape_full_native_slurm.sh
```

The site may require additional `--account` or `--partition` options.
Use an allocated compute node or Slurm batch job, not a login node for
the actual dense processing. `%8` caps simultaneous jobs at eight,
subject to site policies and actual resources. Check `squeue -u $USER`,
`sacct`, and the task logs for scheduling or runtime errors.

The NPZ flags are: 0 finite optimum; 1 optimizer/range warning;
2 unusable background, nonfinite data or fit failure. The seven numeric
fields per channel are `lnL0`, `lnLmax`, signed `mu_hat`, nominal
one-sided `q0`, `max_T/B`, `N_avg` at center, and number of background
segments used. A missing sigma is not an error here: no Fisher
uncertainty is attempted during the dense pass. Do not calculate an
upper limit from the uncalibrated output. Candidate re-analysis,
coverage calibration and physical X/Y response are the next stage.

**Provenance and reruns:** production inputs must not change during the scan.
If a task file already exists it is skipped; do not reuse the same
output directory after changing frequency/grid/physics parameters.
Rename the directory or use a new version for new configurations.

**Resource warning:** nearly five million hypotheses require significant
CPU despite batching. The 200-mass pilot measures the actual throughput
on VIPER and determines whether `--task-bins 10000`, 8-hour walltime,
and array concurrency are appropriate. No claim is made that the full
array has already been submitted or completed.
