#!/usr/bin/env bash
# Generate band plots (SOC/NoSOC, strains 099/100/101) and write a CSV summary.
# Legacy-friendly: no -o, no "vbm" literal; move default band.png outward.

set -euo pipefail

STRAINS="099 100 101"
MODES="soc nosoc"
OUTDIR="sumo_plots/bands"
SUMMARY="sumo_plots/band_summary.csv"
YMIN="-3"
YMAX="3"
FMT="png"

mkdir -p "$OUTDIR"
mkdir -p "sumo_plots"
echo "strain,mode,Eg(eV),type,VBM_k,CBM_k,VBM(eV)" > "$SUMMARY"

get_gap_json() {
python3 - "$1" <<'PY'
import sys,json
from pymatgen.io.vasp.outputs import BSVasprun
vasp = BSVasprun(sys.argv[1], parse_projected_eigen=False)
bs = vasp.get_band_structure(line_mode=True)

gap = bs.get_band_gap()
Eg = float(gap.get("energy", 0.0))
gtype = "direct" if gap.get("direct", False) else "indirect"
vbm = bs.get_vbm()
cbm = bs.get_cbm()

def kp_name(edge):
    kp = edge.get("kpoint")
    if kp is not None and getattr(kp,"label",None):
        return kp.label
    idx = edge.get("kpoint_index")
    return f"k#{idx}" if idx is not None else "unknown"

print(json.dumps({
    "Eg": Eg,
    "type": gtype,
    "vbm_e": float(vbm["energy"]),
    "vbm_k": kp_name(vbm),
    "cbm_k": kp_name(cbm),
}))
PY
}

for s in $STRAINS; do
  for m in $MODES; do
    case "$m" in
      soc)   W="strain_${s}/24_band_soc_HR" ;;
      nosoc) W="strain_${s}/24_band_nosoc_HR" ;;
    esac
    [[ -d "$W" ]] || { echo "[SKIP] $W"; continue; }
    [[ -s "$W/vasprun.xml" && -s "$W/KPOINTS" ]] || { echo "[WARN] missing files in $W"; continue; }

    echo "[BAND] $W → $OUTDIR/band_${s}_${m}.${FMT}"
    ( cd "$W"
      rm -f band.$FMT bandstructure.$FMT BandStructure.$FMT 2>/dev/null || true
      # 레거시: EF=0 + zero-line만. (버전 따라 --zero-energy 숫자 지원 불안정)
      sumo-bandplot -f vasprun.xml -p KPOINTS --zero-line --ymin "$YMIN" --ymax "$YMAX" --format "$FMT"
    )

    # 기본 출력 파일명 탐지 후 이동
    SRC=""
    for cand in "$W/band.$FMT" "$W/bandstructure.$FMT" "$W/BandStructure.$FMT"; do
      [[ -s "$cand" ]] && { SRC="$cand"; break; }
    done
    if [[ -z "$SRC" ]]; then
      # 혹시 다른 이름으로 나왔으면 최근 png 하나 집기
      LATEST=$(ls -t "$W"/*."$FMT" 2>/dev/null | head -n1 || true)
      [[ -n "$LATEST" ]] && SRC="$LATEST"
    fi
    if [[ -n "$SRC" ]]; then
      cp -f "$SRC" "$OUTDIR/band_${s}_${m}.$FMT"
    else
      echo "[WARN] No band image in $W"
    fi

    META=$(get_gap_json "$W/vasprun.xml")
    EG=$(python3 -c "import json; d=json.loads('''$META'''); print(f'{d[\"Eg\"]:.3f}')")
    TYP=$(python3 -c "import json; d=json.loads('''$META'''); print(d['type'])")
    VBMK=$(python3 -c "import json; d=json.loads('''$META'''); print(d['vbm_k'])")
    CBMK=$(python3 -c "import json; d=json.loads('''$META'''); print(d['cbm_k'])")
    VBME=$(python3 -c "import json; d=json.loads('''$META'''); print(f'{d[\"vbm_e\"]:.6f}')")

    echo "strain_${s},${m},${EG},${TYP},${VBMK},${CBMK},${VBME}" >> "$SUMMARY"
  done
done

echo "[DONE] Bands -> $OUTDIR"
echo "[DONE] CSV   -> $SUMMARY"
