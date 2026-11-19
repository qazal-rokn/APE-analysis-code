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
# Ensure the entire requested native frequency interval was actually run.
# Vendor part combiner checks internal continuity but not the final endpoint.
"$PY" - "$DIR" "${APE_FMIN:-65}" "${APE_FMAX:-8000}" <<'PY'
import json, pathlib, sys, h5py, numpy as np
root = pathlib.Path(sys.argv[1])
fmin, fmax = map(float, sys.argv[2:4])
spectra = list(root.glob("APE_*_Q.h5"))
if len(spectra) != 1:
    raise SystemExit(f"Expected exactly one Q spectrum; found {len(spectra)}")
with h5py.File(spectra[0]) as h:
    f = h["frequency"][:]
expected_start = int(np.searchsorted(f, fmin, side="left"))
expected_stop = int(np.searchsorted(f, fmax, side="left"))
manifest_files = sorted((root / "likelihood/parts").glob("Q_*_result.json"))
intervals = []
for path in manifest_files:
    with open(path) as src:
        m = json.load(src)
    indices = m.get("native_grid_indices")
    if m.get("channel") != "Q" or not isinstance(indices, list) or len(indices) != 2:
        raise SystemExit(f"Invalid Q likelihood manifest {path}")
    if not pathlib.Path(str(path).replace(".json",".npy")).is_file():
        raise SystemExit(f"Likelihood data file missing for {path}")
    intervals.append((int(indices[0]), int(indices[1])))
intervals.sort()
if not intervals or intervals[0][0] != expected_start or intervals[-1][1] != expected_stop:
    raise SystemExit(f"Incomplete likelihood: expected [{expected_start}, {expected_stop}), "
                     f"have {intervals[0] if intervals else None} through {intervals[-1] if intervals else None}")
for a, b in zip(intervals, intervals[1:]):
    if a[1] != b[0]:
        raise SystemExit(f"Gap/overlap in likelihood: {a} and {b}")
print(f"Full native-frequency coverage verified: {expected_stop-expected_start} bins")
PY
"$PY" -m ape_post_lpsd.ape_combine_parts \
  --results-dir "$PARTS" --prefix Q --output-npy "$OUT"
[[ -s "$OUT" && -s "${OUT%.npy}.json" ]] || { echo "MLE combiner outputs missing" >&2; exit 3; }
"$PY" -m ape_post_lpsd.ape_finalize_inference \
  --maxlkl-npy "$OUT" \
  --step4-manifest-json "${OUT%.npy}.json" \
  --output-prefix "$PREFIX" --alpha-cl 0.95 --sigma-limit 5
echo "Expected UL: ${PREFIX}_expected_ul.csv"
echo "Candidate screening: ${PREFIX}_screened_candidates.csv"
