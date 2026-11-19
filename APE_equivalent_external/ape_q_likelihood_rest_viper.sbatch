#!/usr/bin/env bash
#SBATCH --job-name=ape-Q-LKL-rest
#SBATCH --partition=small
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --time=12:00:00
#SBATCH --output=ape_Q_likelihood_rest_%A_%a.out
#SBATCH --error=ape_Q_likelihood_rest_%A_%a.err

# APE scheduling-only wrapper: each Slurm task processes several distinct
# original-vendor 1000-bin likelihood chunks, reducing group job-count pressure.
# This does NOT change original likelihood, response, background or outputs.
set -euo pipefail
ROOT="${APE_REPO:-$PWD}"
if [[ -d "$ROOT/APE_equivalent_external" ]]; then ROOT="$ROOT/APE_equivalent_external"; fi
cd "$ROOT"
export PYTHONPATH="$ROOT${PYTHONPATH:+:$PYTHONPATH}"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
PY="${APE_PYTHON:-/ptmp/qarokn/ape-python-venv-20261023/bin/python}"
DIR="${APE_OUTPUT_DIR:-/ptmp/qarokn/APE_Q_FULL_BACKGROUND_20261119}/combined_Q"
Q="$DIR/APE_${APE_T0:-0}_${APE_T1:-15435}_Q.h5"
BKG="$DIR/background_Q.json"
RESP="$DIR/ape_Q_response.npz"
[[ -x "$PY" && -s "$Q" && -s "$BKG" && -s "$RESP" ]] || {
  echo "Missing Python or Q/background/response inputs" >&2; exit 2;
}
[[ -v SLURM_ARRAY_TASK_ID ]] || { echo "Only run as Slurm array" >&2; exit 2; }
START_TASK="${APE_START_TASK:?First remaining chunk index required}"
CHUNKS_PER_TASK="${APE_CHUNKS_PER_TASK:-20}"
TOTAL_CHUNKS="${APE_TOTAL_CHUNKS:?Total chunk count required}"
CHUNK_SIZE="${APE_CHUNK_SIZE:-1000}"
for value in "$START_TASK" "$CHUNKS_PER_TASK" "$TOTAL_CHUNKS" "$CHUNK_SIZE"; do
  [[ "$value" =~ ^[0-9]+$ ]] || { echo "Invalid nonnegative integer: $value" >&2; exit 2; }
done
[[ "$CHUNKS_PER_TASK" -gt 0 && "$CHUNK_SIZE" -gt 0 ]] || exit 2
FIRST=$((START_TASK + SLURM_ARRAY_TASK_ID * CHUNKS_PER_TASK))
END=$((FIRST + CHUNKS_PER_TASK))
if ((END > TOTAL_CHUNKS)); then END="$TOTAL_CHUNKS"; fi
mkdir -p "$DIR/likelihood/parts"
echo "Array index $SLURM_ARRAY_TASK_ID: original-vendor chunks [$FIRST,$END)"
for ((TASK_INDEX=FIRST; TASK_INDEX<END; TASK_INDEX++)); do
  printf -v TASK '%06d' "$TASK_INDEX"
  OUTPUT="$DIR/likelihood/parts/Q_${TASK}_result.npy"
  REPORT="${OUTPUT%.npy}.json"
  if [[ -s "$OUTPUT" && -s "$REPORT" ]]; then
    echo "Already present, not overwriting chunk $TASK_INDEX: $OUTPUT"
    continue
  fi
  if [[ -e "$OUTPUT" || -e "$REPORT" ]]; then
    echo "Partial or incomplete existing result at chunk $TASK_INDEX; stopping, not overwriting." >&2
    exit 3
  fi
  echo "Computing original-vendor likelihood chunk $TASK_INDEX"
  /usr/bin/time -p "$PY" -m configuration_A.run likelihood \
    --spectrum "$Q" --background-json "$BKG" --response-npz "$RESP" \
    --output-npy "$OUTPUT" --ana-fmin "${APE_FMIN:-65}" \
    --ana-fmax "${APE_FMAX:-8000}" \
    --min-log10mu "${APE_MIN_LOG10MU:--8}" --max-log10mu "${APE_MAX_LOG10MU:-8}" \
    --start-index "$((TASK_INDEX * CHUNK_SIZE))" --max-frequencies "$CHUNK_SIZE"
  [[ -s "$OUTPUT" && -s "$REPORT" ]] || { echo "Chunk $TASK_INDEX produced no full result" >&2; exit 3; }
done
echo "All assigned original-vendor likelihood chunks completed."
