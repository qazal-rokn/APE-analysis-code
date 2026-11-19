#!/usr/bin/env bash
#SBATCH --job-name=ape-Q-fastcheck
#SBATCH --partition=small
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:45:00
#SBATCH --output=ape_Q_fastcheck_%j.out
#SBATCH --error=ape_Q_fastcheck_%j.err
# Pilot-validated analytic equivalent of ONE-bin original vendor likelihood.
# Experimental acceleration only; never writes to existing original-job paths.
set -euo pipefail
cd /ptmp/qarokn/APE-analysis-code-configA/APE_equivalent_external
export PYTHONPATH="$PWD"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
PY=/ptmp/qarokn/ape-python-venv-20261023/bin/python
DIR=/ptmp/qarokn/APE_Q_FULL_BACKGROUND_20261119/combined_Q
mkdir -p "$DIR/likelihood/fast_diagnostic"
echo "Starting fast single-bin analytic diagnostic; original vendor jobs untouched."
/usr/bin/time -p "$PY" -u -m configuration_A.fast_q_singlebin_diagnostic \
  --spectrum "$DIR/APE_0_15435_Q.h5" \
  --background-json "$DIR/background_Q.json" \
  --response-npz "$DIR/ape_Q_response.npz" \
  --pilot-npy "$DIR/likelihood/pilot/pilot_000000_result.npy" \
  --output-npy "$DIR/likelihood/fast_diagnostic/Q_full_analytic_diagnostic.npy"
echo "Analytic diagnostic completed"
