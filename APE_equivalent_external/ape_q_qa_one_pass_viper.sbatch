#!/usr/bin/env bash
#SBATCH --job-name=ape-Q-QA
#SBATCH --partition=small
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --time=01:00:00
#SBATCH --output=ape_Q_QA_%j.out
#SBATCH --error=ape_Q_QA_%j.err
set -euo pipefail
cd /ptmp/qarokn/APE-analysis-code-configA/APE_equivalent_external
export PYTHONPATH="$PWD"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
PY=/ptmp/qarokn/ape-python-venv-20261023/bin/python
DIR=/ptmp/qarokn/APE_Q_FULL_BACKGROUND_20261119/combined_Q
mkdir -p "$DIR/qa_one_pass"
echo "One-pass QA for existing 4.2875h PSD/background; no inputs overwritten."
/usr/bin/time -p "$PY" -u -m configuration_A.qa_q_background_all_in_one \
  --spectrum "$DIR/APE_0_15435_Q.h5" \
  --background-json "$DIR/background_Q.json" \
  --output-dir "$DIR/qa_one_pass" --refit
echo "QA completed; outputs in $DIR/qa_one_pass"
