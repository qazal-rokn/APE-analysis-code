#!/usr/bin/env bash
# Submit full 65–8000 Hz likelihood after reviewing a successful 20-bin pilot.
# The numerical likelihood is unchanged vendor calc_max_lkl.
set -euo pipefail
ROOT="${APE_REPO:-$PWD}"
if [[ -d "$ROOT/APE_equivalent_external" ]]; then ROOT="$ROOT/APE_equivalent_external"; fi
cd "$ROOT"
PY="${APE_PYTHON:-/ptmp/qarokn/ape-python-venv-20261023/bin/python}"
DIR="${APE_OUTPUT_DIR:-/ptmp/qarokn/APE_Q_FULL_BACKGROUND_20261119}/combined_Q"
SPECTRUM="$DIR/APE_${APE_T0:-0}_${APE_T1:-15435}_Q.h5"
PILOT="$DIR/likelihood/pilot/pilot_000000_result.npy"
CHUNK="${APE_CHUNK_SIZE:-1000}"
WAVE="${APE_WAVE_SIZE:-250}"
CONCURRENT="${APE_ARRAY_CONCURRENCY:-8}"
[[ -s "$SPECTRUM" && -s "$PILOT" && -s "$DIR/ape_Q_response.npz" ]] || {
  echo "Need Q spectrum, successful pilot and APE response first" >&2; exit 2;
}
[[ "$CHUNK" =~ ^[1-9][0-9]*$ && "$WAVE" =~ ^[1-9][0-9]*$ && "$CONCURRENT" =~ ^[1-9][0-9]*$ ]] || {
  echo "Invalid positive integer CHUNK/WAVE/CONCURRENT" >&2; exit 2;
}
read -r TASKS BINS < <("$PY" - "$SPECTRUM" "$PILOT" "$CHUNK" "${APE_FMIN:-65}" "${APE_FMAX:-8000}" <<'PY'
import sys, math, h5py, numpy as np
spec, pilot, chunk, low, high = sys.argv[1:]
with h5py.File(spec) as h:
    f = h['frequency'][:]
start = int(np.searchsorted(f, float(low), side='left'))
stop = int(np.searchsorted(f, float(high), side='left'))
n = stop-start
p = np.load(pilot, allow_pickle=False)
if p.ndim != 2 or p.shape[1] != 5 or len(p) == 0:
    raise SystemExit('Invalid pilot result')
if not np.any(np.isfinite(p[:, 3]) & np.isfinite(p[:, 4]) & (p[:,4] > 0)):
    raise SystemExit('Pilot has no valid mu/sigma fits. Do not scale up blindly.')
if n <= 0:
    raise SystemExit('No frequency bins in selected band')
print(math.ceil(n/int(chunk)), n)
PY
)
echo "Requested $BINS frequency bins in $TASKS likelihood chunks."
echo "Each original-vendor optimizer processes one bin; this may take much longer than today."
for ((first=0; first<TASKS; first+=WAVE)); do
  count=$((TASKS-first))
  if ((count>WAVE)); then count="$WAVE"; fi
  last=$((count-1))
  echo "Submitting task indices $first through $((first+last))"
  sbatch --parsable --array="0-${last}%${CONCURRENT}" \
    --export="ALL,APE_REPO=$ROOT,APE_LKL_MODE=full,APE_CHUNK_SIZE=$CHUNK,APE_START_TASK=$first" \
    ape_q_likelihood_viper.sbatch
done
echo "After ALL array tasks finish successfully, submit ape_q_finalize_viper.sbatch."
