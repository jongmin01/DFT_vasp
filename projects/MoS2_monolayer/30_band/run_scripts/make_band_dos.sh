#!/bin/bash
# make_band_dos_v2.sh
# SOC/NOSOC 디렉토리가 뒤바뀌어 읽히는 문제를 완전히 차단
# 실행 로그를 명확히 찍어서 사람이 보면서도 검증 가능

set -euo pipefail

ROOT="$(pwd)"
OUT_BANDS="$ROOT/plots/bands"
OUT_DOS="$ROOT/plots/dos"
OUT_SUM="$ROOT/plots/summary_band.csv"

EMIN="-3.0"
EMAX="3.0"
DOS_XMIN="-5.0"
DOS_XMAX="5.0"

mkdir -p "$OUT_BANDS" "$OUT_DOS"

xml_ok () {
  local f="$1"
  [[ -s "$f" ]] || return 1
  tail -n 1 "$f" | grep -q '</modeling>' && return 0 || return 1
}

py_gap_meta () {
python3 - "$1" <<'PY'
import sys, json
from pymatgen.io.vasp.outputs import BSVasprun
vasp_xml = sys.argv[1]
v = BSVasprun(vasp_xml, parse_projected_eigen=False)
bs = v.get_band_structure(line_mode=True)
gapinfo = bs.get_band_gap()
Eg = float(gapinfo.get("energy",0.0))
gtype = "direct" if gapinfo.get("direct", False) else "indirect"
def lbl(edge):
    if hasattr(edge.get("kpoint",""), "label") and edge["kpoint"].label:
        return edge["kpoint"].label
    if "kpoint_index" in edge and edge["kpoint_index"] is not None:
        return f"k#{edge['kpoint_index']}"
    return "unknown"
vbm_e = float(bs.get_vbm()["energy"])
vbm_k = lbl(bs.get_vbm())
cbm_k = lbl(bs.get_cbm())
print(json.dumps({"Eg":Eg, "type":gtype, "vbm_e":vbm_e, "vbm_k":vbm_k, "cbm_k":cbm_k}))
PY
}

py_get_ef () {
python3 - "$1" <<'PY'
import sys
from pymatgen.io.vasp.outputs import Vasprun
try:
    v = Vasprun(sys.argv[1], parse_potcar_file=False)
    print(f"{v.efermi:.6f}")
except:
    print("")
PY
}

echo "strain,mode,Eg(eV),type,VBM_k,CBM_k,VBM_energy(eV),band_png,dos_png" > "$OUT_SUM"

found_any=false

for S in $(ls -d strain_* 2>/dev/null | sort); do
  found_any=true

  # ======== SOC 먼저 ========
  MODE="soc"
  BAND_DIR="$S/24_band_soc_HR"
  DOS_DIR="$S/12_scf_soc_hr"
  TAG="${S#strain_}_soc"

  BAND_XML="$BAND_DIR/vasprun.xml"
  DOS_XML="$DOS_DIR/vasprun.xml"

  echo
  echo "======================"
  echo "[SOC] strain = $S"
  echo "[DEBUG] BAND_DIR = $BAND_DIR"
  echo "[DEBUG] DOS_DIR  = $DOS_DIR"
  echo "======================"

  band_png_rel=""
  dos_png_rel=""

  if xml_ok "$BAND_XML"; then
    echo "[SOC-BAND] plotting $BAND_XML"
    python3 "$ROOT/plot_band_from_vasprun.py" -f "$BAND_XML" \
      -o "$OUT_BANDS/band_${TAG}.png" \
      --emin "$EMIN" --emax "$EMAX" \
      --mark-edges --label-edges \
    && band_png_rel="plots/bands/band_${TAG}.png"
  else
    echo "[WARN] SOC BAND XML unusable: $BAND_XML"
  fi

  echo "[SOC-DOS] plotting $DOS_XML"
  if xml_ok "$DOS_XML"; then
    EF="$(py_get_ef "$DOS_XML")"
    echo "[DEBUG] SOC EF = $EF"
    cd "$DOS_DIR"
    if [[ -n "$EF" ]]; then
      sumo-dosplot -f vasprun.xml \
        --zero-energy "$EF" \
        --xmin "$DOS_XMIN" --xmax "$DOS_XMAX" \
        --format png -o "$OUT_DOS/dos_${TAG}.png"
    else
      sumo-dosplot -f vasprun.xml \
        --no-shift --zero-line \
        --xmin "$DOS_XMIN" --xmax "$DOS_XMAX" \
        --format png -o "$OUT_DOS/dos_${TAG}.png"
    fi
    dos_png_rel="plots/dos/dos_${TAG}.png"
    cd "$ROOT"
  else
    echo "[WARN] SOC DOS XML unusable: $DOS_XML"
  fi

  # gap summary
  EG=""; GTYPE=""; VBMK=""; CBMK=""; VBME=""
  if xml_ok "$BAND_XML"; then
    META_JSON="$(py_gap_meta "$BAND_XML" 2>/dev/null || true)"
    if [[ -n "$META_JSON" ]]; then
      EG="$(python3 -c "import json; d=json.loads('''$META_JSON'''); print(f'{d[\"Eg\"]:.3f}')" 2>/dev/null || echo "")"
      GTYPE="$(python3 -c "import json; d=json.loads('''$META_JSON'''); print(d['type'])" 2>/dev/null || echo "")"
      VBMK="$(python3 -c "import json; d=json.loads('''$META_JSON'''); print(d['vbm_k'])" 2>/dev/null || echo "")"
      CBMK="$(python3 -c "import json; d=json.loads('''$META_JSON'''); print(d['cbm_k'])" 2>/dev/null || echo "")"
      VBME="$(python3 -c "import json; d=json.loads('''$META_JSON'''); print(f'{d[\"vbm_e\"]:.6f}')" 2>/dev/null || echo "")"
    fi
  fi

  echo "${S},SOC,${EG},${GTYPE},${VBMK},${CBMK},${VBME},${band_png_rel},${dos_png_rel}" >> "$OUT_SUM"



  # ======== NOSOC ========
  MODE="nosoc"
  BAND_DIR="$S/24_band_nosoc_HR"
  DOS_DIR="$S/12_scf_nosoc_hr"
  TAG="${S#strain_}_nosoc"

  BAND_XML="$BAND_DIR/vasprun.xml"
  DOS_XML="$DOS_DIR/vasprun.xml"

  echo
  echo "======================"
  echo "[NOSOC] strain = $S"
  echo "[DEBUG] BAND_DIR = $BAND_DIR"
  echo "[DEBUG] DOS_DIR  = $DOS_DIR"
  echo "======================"

  band_png_rel=""
  dos_png_rel=""

  if xml_ok "$BAND_XML"; then
    echo "[NOSOC-BAND] plotting $BAND_XML"
    python3 "$ROOT/plot_band_from_vasprun.py" -f "$BAND_XML" \
      -o "$OUT_BANDS/band_${TAG}.png" \
      --emin "$EMIN" --emax "$EMAX" \
      --mark-edges --label-edges \
    && band_png_rel="plots/bands/band_${TAG}.png"
  else
    echo "[WARN] NOSOC BAND XML unusable: $BAND_XML"
  fi

  echo "[NOSOC-DOS] plotting $DOS_XML"
  if xml_ok "$DOS_XML"; then
    EF="$(py_get_ef "$DOS_XML")"
    echo "[DEBUG] NOSOC EF = $EF"
    cd "$DOS_DIR"
    if [[ -n "$EF" ]]; then
      sumo-dosplot -f vasprun.xml \
        --zero-energy "$EF" \
        --xmin "$DOS_XMIN" --xmax "$DOS_XMAX" \
        --format png -o "$OUT_DOS/dos_${TAG}.png"
    else
      sumo-dosplot -f vasprun.xml \
        --no-shift --zero-line \
        --xmin "$DOS_XMIN" --xmax "$DOS_XMAX" \
        --format png -o "$OUT_DOS/dos_${TAG}.png"
    fi
    dos_png_rel="plots/dos/dos_${TAG}.png"
    cd "$ROOT"
  else
    echo "[WARN] NOSOC DOS XML unusable: $DOS_XML"
  fi

  # gap summary
  EG=""; GTYPE=""; VBMK=""; CBMK=""; VBME=""
  if xml_ok "$BAND_XML"; then
    META_JSON="$(py_gap_meta "$BAND_XML" 2>/dev/null || true)"
    if [[ -n "$META_JSON" ]]; then
      EG="$(python3 -c "import json; d=json.loads('''$META_JSON'''); print(f'{d[\"Eg\"]:.3f}')" 2>/dev/null || echo "")"
      GTYPE="$(python3 -c "import json; d=json.loads('''$META_JSON'''); print(d['type'])" 2>/dev/null || echo "")"
      VBMK="$(python3 -c "import json; d=json.loads('''$META_JSON'''); print(d['vbm_k'])" 2>/dev/null || echo "")"
      CBMK="$(python3 -c "import json; d=json.loads('''$META_JSON'''); print(d['cbm_k'])" 2>/dev/null || echo "")"
      VBME="$(python3 -c "import json; d=json.loads('''$META_JSON'''); print(f'{d[\"vbm_e\"]:.6f}')" 2>/dev/null || echo "")"
    fi
  fi

  echo "${S},NOSOC,${EG},${GTYPE},${VBMK},${CBMK},${VBME},${band_png_rel},${dos_png_rel}" >> "$OUT_SUM"

done

if ! $found_any; then
  echo "[ERROR] strain_* 디렉토리 없음. MoS2 루트에서 실행하세요."
  exit 1
fi

echo
echo "[DONE] Band images : $OUT_BANDS"
echo "[DONE] DOS images  : $OUT_DOS"
echo "[DONE] Summary CSV : $OUT_SUM"

