#!/usr/bin/env bash
#SBATCH --job-name=ape-Q-likelihood
#SBATCH --partition=small
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --time=12:00:00
#SBATCH --output=ape_Q_likelihood_%A_%a.out
#SBATCH --error=ape_Q_likelihood_%A_%a.err

# One unchanged-vendor optimizer per frequency, driven by Configuration A.
# Pilot first: sbatch --export=ALL,APE_LKL_MODE=pilot ape_q_likelihood_viper.sbatch
# Full: sbatch --array=0-N%W --export=ALL,APE_LKL_MODE=full,APE_CHUNK_SIZE=1000 ...
set -euo pipefail
ROOT="${APE_REPO:-$PWD}"
if [[ -d "$ROOT/APE_equivalent_external" ]]; then ROOT="$ROOT/APE_equivalent_external"; fi
cd "$ROOT"
export PYTHONPATH="$ROOT${PYTHONPATH:+:$PYTHONPATH}"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
PY="${APE_PYTHON:-/ptmp/qarokn/ape-python-venv-20261023/bin/python}"
[[ -x "$PY" ]] || { echo "APE_PYTHON not found: $PY" >&2; exit 2; }
"$PY" -c 'import numpy, scipy, h5py, pandas, findiff, scalar_dark_matter_lpsd'
DIR="${APE_OUTPUT_DIR:-/ptmp/qarokn/APE_Q_FULL_BACKGROUND_20261119}/combined_Q"
Q="$DIR/APE_${APE_T0:-0}_${APE_T1:-15435}_Q.h5"
BKG="$DIR/background_Q.json"
RESP="$DIR/ape_Q_response.npz"
[[ -s "$Q" && -s "$BKG" ]] || { echo "Q background/spectrum not ready: $Q / $BKG" >&2; exit 2; }
MODE="${APE_LKL_MODE:-pilot}"
FMIN="${APE_FMIN:-65}"; FMAX="${APE_FMAX:-8000}"
MINMU="${APE_MIN_LOG10MU:--8}"; MAXMU="${APE_MAX_LOG10MU:-8}"
if [[ "$MODE" == "pilot" ]]; then
  if [[ ! -s "$RESP" ]]; then
    "$PY" -m configuration_A.run response --spectrum "$Q" --output-npz "$RESP" \
      --finesse-s "${APE_FINESSE_S:-1500}" --tau-s 2.5e-9 --tau1-s 1.0e-8
  fi
  mkdir -p "$DIR/likelihood/pilot"
  OUTPUT="$DIR/likelihood/pilot/pilot_000000_result.npy"
  START=0
  COUNT="${APE_PILOT_BINS:-20}"
elif [[ "$MODE" == "full" ]]; then
  [[ -s "$RESP" ]] || { echo "Generate response by successful pilot first" >&2; exit 2; }
  [[ -v SLURM_ARRAY_TASK_ID ]] || { echo "Full mode requires a Slurm array index" >&2; exit 2; }
  CHUNK="${APE_CHUNK_SIZE:-1000}"
  TASK_INDEX=$(( ${APE_START_TASK:-0} + SLURM_ARRAY_TASK_ID ))
  START=$((TASK_INDEX * CHUNK))
  COUNT="$CHUNK"
  mkdir -p "$DIR/likelihood/parts"
  printf -v TASK '%06d' "$TASK_INDEX"
  OUTPUT="$DIR/likelihood/parts/Q_${TASK}_result.npy"
else
  echo "APE_LKL_MODE must be pilot or full" >&2; exit 2
fi
[[ ! -e "$OUTPUT" ]] || { echo "Refusing existing output $OUTPUT" >&2; exit 2; }
echo "mode=$MODE start_index=$START count=$COUNT bounds=$MINMU,$MAXMU"
/usr/bin/time -p "$PY" -m configuration_A.run likelihood \
  --spectrum "$Q" --background-json "$BKG" --response-npz "$RESP" \
  --output-npy "$OUTPUT" --ana-fmin "$FMIN" --ana-fmax "$FMAX" \
  --min-log10mu "$MINMU" --max-log10mu "$MAXMU" \
  --start-index "$START" --max-frequencies "$COUNT"
echo "Likelihood segment saved: $OUTPUT"
