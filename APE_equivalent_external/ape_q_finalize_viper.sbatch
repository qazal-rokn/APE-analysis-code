#!/usr/bin/env bash
#SBATCH --job-name=ape-Q-finalize
#SBATCH --partition=small
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=ape_Q_finalize_%j.out
#SBATCH --error=ape_Q_finalize_%j.err

# Run ONLY after every full likelihood array task is COMPLETED successfully.
# Original vendor numerical combiner and original expected-UL method.
set -euo pipefail
ROOT="${APE_REPO:-$PWD}"
if [[ -d "$ROOT/APE_equivalent_external" ]]; then ROOT="$ROOT/APE_equivalent_external"; fi
cd "$ROOT"
export PYTHONPATH="$ROOT${PYTHONPATH:+:$PYTHONPATH}"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
PY="${APE_PYTHON:-/ptmp/qarokn/ape-python-venv-20261023/bin/python}"
[[ -x "$PY" ]] || { echo "APE_PYTHON missing" >&2; exit 2; }
DIR="${APE_OUTPUT_DIR:-/ptmp/qarokn/APE_Q_FULL_BACKGROUND_20261119}/combined_Q"
PARTS="$DIR/likelihood/parts"
OUT="$DIR/likelihood/likelihood_Q_full.npy"
PREFIX="$DIR/likelihood/ape_Q"
[[ -d "$PARTS" ]] || { echo "Likelihood part directory not present" >&2; exit 2; }
"$PY" -m ape_post_lpsd.ape_combine_parts \
  --results-dir "$PARTS" --prefix Q --output-npy "$OUT"
[[ -s "$OUT" && -s "${OUT%.npy}.json" ]] || { echo "MLE combiner outputs missing" >&2; exit 3; }
"$PY" -m ape_post_lpsd.ape_finalize_inference \
  --maxlkl-npy "$OUT" \
  --step4-manifest-json "${OUT%.npy}.json" \
  --output-prefix "$PREFIX" --alpha-cl 0.95 --sigma-limit 5
echo "Expected UL: ${PREFIX}_expected_ul.csv"
echo "Candidate screening: ${PREFIX}_screened_candidates.csv"
