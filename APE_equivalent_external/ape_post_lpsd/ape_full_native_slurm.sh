#!/bin/bash
#SBATCH --job-name=APE_4h_native_LL
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=08:00:00
# Submit with --array and --output; preserve site-specific partition/account.

set -euo pipefail
: "${ROOT:?Set ROOT with sbatch --export}"
: "${OUT:?Set OUT with sbatch --export}"
: "${SLURM_ARRAY_TASK_ID:?Must submit as a Slurm array}"

cd /ptmp/qarokn/APE-analysis-code/APE_equivalent_external
export PYTHONPATH="$PWD"
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

exec /ptmp/qarokn/ape-python-venv-20261023/bin/python \
  -m ape_post_lpsd.ape_full_native_likelihood_scan \
  --root "$ROOT" \
  --outdir "$OUT" \
  --task-bins 10000 \
  --task-id "$SLURM_ARRAY_TASK_ID" \
  --fmin 65.1 --fmax 8000
