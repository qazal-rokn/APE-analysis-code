#!/usr/bin/env bash
#SBATCH --job-name=ape-Q-LKL-8CPU
#SBATCH --partition=small
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=12:00:00
#SBATCH --output=ape_Q_likelihood_8cpu_%A_%a.out
#SBATCH --error=ape_Q_likelihood_8cpu_%A_%a.err

# Use 8 allocated CPUs for 8 independent unchanged-vendor likelihood runs.
# Only scheduling/parallelization differs from the successful 1-core pilot.
# Original response, PSD, background JSON and optimizer stay untouched.
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

[[ "${SLURM_CPUS_PER_TASK:-0}" -eq 8 ]] || { echo "Require exactly 8 allocated CPUs" >&2; exit 2; }
[[ -v SLURM_ARRAY_TASK_ID ]] || { echo "Require a Slurm array index" >&2; exit 2; }
[[ -x "$PY" && -s "$Q" && -s "$BKG" && -s "$RESP" ]] || {
  echo "Missing Python, Q spectrum, fitted background or response" >&2; exit 2;
}
START_TASK="${APE_START_TASK:?Specify first global 1000-bin chunk index}"
TOTAL_CHUNKS="${APE_TOTAL_CHUNKS:?Specify total chunks}"
CHUNKS_PER_TASK="${APE_CHUNKS_PER_TASK:-50}"
CHUNK_SIZE="${APE_CHUNK_SIZE:-1000}"
for value in "$START_TASK" "$TOTAL_CHUNKS" "$CHUNKS_PER_TASK" "$CHUNK_SIZE"; do
  [[ "$value" =~ ^[0-9]+$ ]] || { echo "Not a nonnegative integer: $value" >&2; exit 2; }
done
[[ "$CHUNKS_PER_TASK" -gt 0 && "$CHUNK_SIZE" -gt 0 ]] || exit 2

FIRST=$((START_TASK + SLURM_ARRAY_TASK_ID * CHUNKS_PER_TASK))
LAST=$((FIRST + CHUNKS_PER_TASK))
if ((LAST > TOTAL_CHUNKS)); then LAST="$TOTAL_CHUNKS"; fi
echo "Executing original-vendor likelihood chunks [$FIRST, $LAST), 8 concurrent one-core workers"
mkdir -p "$DIR/likelihood/parts" "$DIR/likelihood/worker_logs"

run_one() {
  local k="$1" task output report log
  printf -v task '%06d' "$k"
  output="$DIR/likelihood/parts/Q_${task}_result.npy"
  report="${output%.npy}.json"
  log="$DIR/likelihood/worker_logs/Q_${task}_job_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.log"

  if [[ -s "$output" && -s "$report" ]]; then
    echo "Already finished; keeping existing chunk $k"
    return 0
  fi
  if [[ -e "$output" || -e "$report" || -e "${output%.npy}.partial.npy" ]]; then
    echo "Incomplete existing result for chunk $k: refusing overwrite" >&2
    return 3
  fi
  echo "Starting chunk $k, writing $log"
  if ! /usr/bin/time -p "$PY" -u -m configuration_A.run likelihood \
    --spectrum "$Q" --background-json "$BKG" --response-npz "$RESP" \
    --output-npy "$output" --ana-fmin "${APE_FMIN:-65}" --ana-fmax "${APE_FMAX:-8000}" \
    --min-log10mu "${APE_MIN_LOG10MU:--8}" --max-log10mu "${APE_MAX_LOG10MU:-8}" \
    --start-index "$((k * CHUNK_SIZE))" --max-frequencies "$CHUNK_SIZE" > "$log" 2>&1
  then
    echo "Vendor likelihood failed in chunk $k (last 15 lines follow)" >&2
    tail -15 "$log" >&2
    return 3
  fi
  [[ -s "$output" && -s "$report" ]] || {
    echo "Chunk $k has incomplete output after vendor run" >&2
    return 3
  }
  echo "Completed chunk $k"
}

pids=()
failed=0
for ((k=FIRST;k<LAST;k++)); do
  run_one "$k" &
  pids+=("$!")
  if (("${#pids[@]}" == 8)); then
    for pid in "${pids[@]}"; do
      wait "$pid" || failed=1
    done
    if ((failed)); then
      echo "At least one worker failed, aborting next group; successful results preserved" >&2
      exit 3
    fi
    pids=()
  fi
done
for pid in "${pids[@]}"; do
  wait "$pid" || failed=1
done
if ((failed)); then
  echo "At least one worker failed; completed results preserved" >&2
  exit 3
fi
echo "All assigned chunks successfully completed"
