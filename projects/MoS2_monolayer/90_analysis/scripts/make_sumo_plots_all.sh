#!/bin/bash
# make_sumo_plots_all.sh (Carbon-compatible version)
# Generate band plots with sumo for all strain_* / {SOC,NoSOC} and summarize gaps.

set -euo pipefail

# -------- Settings --------
YMIN="-3"
YMAX="3"
IMG_FMT="png"
OUTDIR="sumo_plots"
SUMMARY="sumo_band_summary.csv"
# ---------------------------

mkdir -p "$OUTDIR"
echo "strain,mode,Eg(eV),type,VBM_k,CBM_k,zero_ref" > "$OUTDIR/$SUMMARY"

# Python helper to get Eg, type, and VBM energy
py_eval_vbm_gap () {
python3 - "$1" <<'PY'
import sys, json
from pymatgen.io.vasp.outputs import BSVasprun

vasp_xml = sys.argv[1]
v = BSVasprun(vasp_xml, parse_projected_eigen=False)
bs = v.get_band_structure(line_mode=True)
vbm_e = bs.get_vbm()["energy"]
gapinfo = bs.get_band_gap()
Eg = float(gapinfo.get("energy", 0.0))
gtype = "direct" if gapinfo.get("direct", False) else "indirect"

def label_or_index(edge):
    if "kpoint" in edge and hasattr(edge["kpoint"], "label") and edge["kpoint"].label:
        return edge["kpoint"].label
    if "kpoint_index" in edge and edge["kpoint_index"] is not None:
        return f"k#{edge['kpoint_index']}"
    return "unknown"

vbm_k = label_or_index(bs.get_vbm())
cbm_k = label_or_index(bs.get_cbm())

print(json.dumps({
    "vbm": vbm_e,
    "Eg": Eg,
    "type": gtype,
    "vbm_k": vbm_k,
    "cbm_k": cbm_k
}))
PY
}

# Loop over strain directories (mapfile replacement)
for S in $(ls -d strain_* 2>/dev/null | sort); do
  for MODE in "soc" "nosoc"; do
    case "$MODE" in
      soc)   WDIR="$S/24_band_soc_HR" ;;
      nosoc) WDIR="$S/24_band_nosoc_HR" ;;
    esac

    if [[ ! -d "$WDIR" ]]; then
      echo "[SKIP] $WDIR not found"; continue
    fi

    VXML="$WDIR/vasprun.xml"
    KPTS="$WDIR/KPOINTS"
    if [[ ! -s "$VXML" ]]; then
      echo "[WARN] $VXML missing → skip $S/$MODE"; continue
    fi
    if [[ ! -s "$KPTS" ]]; then
      echo "[WARN] $KPTS missing → skip $S/$MODE"; continue
    fi

    META_JSON=$(py_eval_vbm_gap "$VXML")
    VBM=$(python3 -c "import json; d=json.loads('''$META_JSON'''); print(f'{d[\"vbm\"]:.6f}')")
    EG=$(python3 -c "import json; d=json.loads('''$META_JSON'''); print(f'{d[\"Eg\"]:.3f}')")
    GTYPE=$(python3 -c "import json; d=json.loads('''$META_JSON'''); print(d['type'])")
    VBMK=$(python3 -c "import json; d=json.loads('''$META_JSON'''); print(d['vbm_k'])")
    CBMK=$(python3 -c "import json; d=json.loads('''$META_JSON'''); print(d['cbm_k'])")

    echo "[RUN] sumo-bandplot for $S ($MODE) with zero=$VBM"
    (
      cd "$WDIR"
      rm -f band.$IMG_FMT bandstructure.$IMG_FMT BandStructure.$IMG_FMT 2>/dev/null || true

      sumo-bandplot -f vasprun.xml -p KPOINTS \
        --zero-energy "$VBM" \
        --ymin "$YMIN" --ymax "$YMAX" \
        --format "$IMG_FMT"
    )

    SRC_IMG=""
    for cand in "$WDIR/band.$IMG_FMT" "$WDIR/bandstructure.$IMG_FMT" "$WDIR/BandStructure.$IMG_FMT"; do
      if [[ -s "$cand" ]]; then SRC_IMG="$cand"; break; fi
    done
    if [[ -z "$SRC_IMG" ]]; then
      LATEST=$(ls -t "$WDIR"/*."$IMG_FMT" 2>/dev/null | head -n1 || true)
      if [[ -n "$LATEST" ]]; then SRC_IMG="$LATEST"; fi
    fi

    if [[ -z "$SRC_IMG" ]]; then
      echo "[WARN] no sumo output image found in $WDIR"
    else
      OUT_IMG="$OUTDIR/${S}_${MODE}_band.$IMG_FMT"
      cp -f "$SRC_IMG" "$OUT_IMG"
      echo "[OK] -> $OUT_IMG"
    fi

    echo "${S},${MODE},${EG},${GTYPE},${VBMK},${CBMK},VBM(${VBM} eV)" >> "$OUTDIR/$SUMMARY"
  done
done

echo
echo "[DONE] Plots saved under: $OUTDIR/"
echo "[DONE] Summary CSV:       $OUTDIR/$SUMMARY"
