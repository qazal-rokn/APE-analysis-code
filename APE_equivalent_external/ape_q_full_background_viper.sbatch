#!/usr/bin/env bash
#SBATCH --job-name=ape-Q-full-bkg
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=23:50:00
#SBATCH --output=ape_Q_full_%j.out
#SBATCH --error=ape_Q_full_%j.err

# Configuration A, unchanged original vendor background fitter.
# Reuses existing 4-column X/Y C-LPSD outputs; sums linear PSDs.
set -euo pipefail
ROOT="${APE_REPO:-$PWD}"
if [[ -d "$ROOT/APE_equivalent_external" ]]; then ROOT="$ROOT/APE_equivalent_external"; fi
[[ -f "$ROOT/configuration_A/run.py" ]] || { echo "APE_REPO does not point to APE_equivalent_external: $ROOT" >&2; exit 2; }
cd "$ROOT"
export PYTHONPATH="$ROOT${PYTHONPATH:+:$PYTHONPATH}"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1

PY="${APE_PYTHON:-$(command -v python3 || command -v python)}"
[[ -x "$PY" ]] || { echo "Set APE_PYTHON to a working Python on Viper" >&2; exit 2; }
"$PY" -c 'import numpy, scipy, h5py, pandas, scalar_dark_matter_lpsd' || {
  echo "Python environment missing project dependencies. Set APE_PYTHON." >&2; exit 2;
}

X="${APE_X_PSD:-/ptmp/qarokn/APE_C_LPSD_pipeline/results/rho_x_combined.lpsd}"
Y="${APE_Y_PSD:-/ptmp/qarokn/APE_C_LPSD_pipeline/results/rho_y_combined.lpsd}"
OUT="${APE_OUTPUT_DIR:-/ptmp/qarokn/APE_Q_FULL_BACKGROUND_20261119}"
T0="${APE_T0:-0}"; T1="${APE_T1:-15435}"
FS="${APE_FS:-20000}"
FMIN="${APE_FMIN:-65}"; FMAX="${APE_FMAX:-8000}"
[[ -s "$X" && -s "$Y" ]] || { echo "Missing X/Y input files; set APE_X_PSD and APE_Y_PSD" >&2; exit 2; }
HANDOFF="$OUT/input_handoff"
QDIR="$OUT/combined_Q"
mkdir -p "$HANDOFF" "$QDIR"
XH="$HANDOFF/APE_${T0}_${T1}_X.h5"
YH="$HANDOFF/APE_${T0}_${T1}_Y.h5"
QH="$QDIR/APE_${T0}_${T1}_Q.h5"
BKG="$QDIR/background_Q.json"

if [[ -s "$XH" && -s "$YH" ]]; then
  echo "Using staged X/Y HDF5 inputs"
elif [[ -e "$XH" || -e "$YH" ]]; then
  echo "Incomplete HDF5 handoff; refusing to overwrite" >&2; exit 3
else
  "$PY" -m ape_post_lpsd.prepare_vendor_lpsd \
    --x-psd "$X" --y-psd "$Y" --output-dir "$HANDOFF" \
    --fs "$FS" --t0 "$T0" --t1 "$T1" --expected-bins 5035956
fi

if [[ -s "$BKG" ]]; then
  echo "Background exists, not overwriting: $BKG"
elif [[ -s "$QH" ]]; then
  "$PY" -m scalar_dark_matter_lpsd.lvk.fit_background \
    --data-path "$QH" --output-path "$BKG" \
    --ana-fmin "$FMIN" --ana-fmax "$FMAX" \
    --n-processes "${SLURM_CPUS_PER_TASK:-8}" --segment-size 1000
else
  "$PY" -m configuration_A.run background \
    --x-h5 "$XH" --y-h5 "$YH" --output-dir "$QDIR" \
    --ana-fmin "$FMIN" --ana-fmax "$FMAX" \
    --n-processes "${SLURM_CPUS_PER_TASK:-8}" --segment-size 1000
fi
[[ -s "$QH" && -s "$BKG" ]] || { echo "Missing expected outputs" >&2; exit 4; }
echo "Completed Q spectrum: $QH"
echo "Completed vendor background fit: $BKG"
