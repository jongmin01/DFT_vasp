#!/usr/bin/env bash
set -euo pipefail
src="../12_scf_nosoc_hr"
for f in CHGCAR WAVECAR POTCAR; do
  if [[ ! -s "$src/$f" ]]; then
    echo "[ERR] $src/$f is missing or empty" >&2
    exit 1
  fi
done
rsync -av "$src/CHGCAR" "$src/WAVECAR" "$src/POTCAR" ./
ls -lh CHGCAR WAVECAR POTCAR
