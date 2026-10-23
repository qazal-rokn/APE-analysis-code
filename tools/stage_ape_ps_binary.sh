#!/usr/bin/env bash
# Stage an already validated C PS source/binary snapshot into a fresh isolated
# directory. Never touches original tarball or vendor source.
set -euo pipefail
if [[ $# -ne 2 ]]; then
  echo "Usage: bash tools/stage_ape_ps_binary.sh ARCHIVE.tar.gz NEW_DEST_DIR" >&2
  exit 2
fi
archive="$1"
destination="$2"
expected="0297b50652b06be90257e7471dc2cfce84cd15ed45add3129fdd34e0c2d6d48e"
[[ -f "$archive" ]] || { echo "Archive missing" >&2; exit 2; }
[[ ! -e "$destination" ]] || { echo "Destination already exists" >&2; exit 2; }
actual=$(shasum -a 256 "$archive" | awk '{print $1}')
[[ "$actual" == "$expected" ]] || { echo "Archive SHA-256 MISMATCH" >&2; exit 1; }
# Reject absolute or parent-traversal archive paths before extraction.
if tar -tzf "$archive" | awk '/^\// || /(^|\/)\.\.($|\/)/ {bad=1} END {exit !bad}'; then
  echo "Unsafe archive path" >&2
  exit 1
fi
mkdir "$destination"
tar -xzf "$archive" -C "$destination"
binary="$destination/build/lpsd-exec"
if [[ ! -f "$binary" ]]; then
  echo "Extracted build/lpsd-exec not found; inspect archive layout" >&2
  exit 1
fi
chmod u+x "$binary"
echo "VERIFIED archive SHA-256: $actual"
shasum -a 256 "$binary"
file "$binary"
echo "Use extracted executable ONLY on a compatible machine (e.g. Linux VIPER)."
