#!/usr/bin/env bash
# Generate DOS plots (SOC/NoSOC, strains 099/100/101). Legacy-friendly.

set -euo pipefail

STRAINS="099 100 101"
MODES="soc nosoc"
OUTDIR="sumo_plots/dos"
FMT="png"
# 필요시 축 조정 (버전에 따라 --ymin/--ymax 동작이 다를 수 있음)
YMIN=""
YMAX=""

mkdir -p "$OUTDIR"

for s in $STRAINS; do
  for m in $MODES; do
    case "$m" in
      soc)   W="strain_${s}/21_dos_soc_HR" ;;
      nosoc) W="strain_${s}/21_dos_nosoc_HR" ;;
    esac
    [[ -d "$W" ]] || { echo "[SKIP] $W"; continue; }
    [[ -s "$W/vasprun.xml" ]] || { echo "[WARN] missing vasprun.xml in $W"; continue; }

    echo "[DOS] $W → $OUTDIR/dos_${s}_${m}.${FMT}"
    ( cd "$W"
      rm -f dos.$FMT 2>/dev/null || true
      # 레거시: 기본 파일명 사용, EF=0 기준
      if [[ -n "$YMIN" && -n "$YMAX" ]]; then
        sumo-dosplot -f vasprun.xml --format "$FMT" --ymin "$YMIN" --ymax "$YMAX"
      else
        sumo-dosplot -f vasprun.xml --format "$FMT"
      fi
    )
    [[ -s "$W/dos.$FMT" ]] && cp -f "$W/dos.$FMT" "$OUTDIR/dos_${s}_${m}.$FMT" || echo "[WARN] No dos image in $W"
  done
done

echo "[DONE] DOS -> $OUTDIR"
