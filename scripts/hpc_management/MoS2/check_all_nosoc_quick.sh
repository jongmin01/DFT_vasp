#!/usr/bin/env bash
# NoSOC summary (EIGENVAL parser with E-fermi from OUTCAR)

set -euo pipefail
shopt -s nullglob

echo "strain, step, converged, nkpts, toten(eV), bandgap(eV), type"

strains=$(ls -d strain_* 2>/dev/null || true)
[ -n "$strains" ] || { echo "[ERR] no strain_* under $(pwd)" >&2; exit 0; }

get_nkpts() {
  awk '/NKPTS/ {for (i=1;i<=NF;i++) if ($i=="NKPTS"){print $(i+2); exit}}' "$1" 2>/dev/null
}

for s in $strains; do
  # step1
  d="$s/11_scf_nosoc_warm"
  if [ -d "$d" ]; then
    conv=$(grep -q reached "$d/OUTCAR" 2>/dev/null && echo YES || echo NO)
    nk=$(get_nkpts "$d/OUTCAR"); nk=${nk:-NA}
    e=$(grep -m1 "free  energy   TOTEN" "$d/OUTCAR" 2>/dev/null | awk '{print $(NF-1)}'); e=${e:-NA}
    echo "$s, step1, $conv, $nk, $e, -, -"
  fi

  # step2
  d="$s/12_scf_nosoc_hr"
  if [ -d "$d" ]; then
    conv=$(grep -q reached "$d/OUTCAR" 2>/dev/null && echo YES || echo NO)
    nk=$(get_nkpts "$d/OUTCAR"); nk=${nk:-NA}
    e=$(grep -m1 "free  energy   TOTEN" "$d/OUTCAR" 2>/dev/null | awk '{print $(NF-1)}'); e=${e:-NA}
    echo "$s, step2, $conv, $nk, $e, -, -"
  fi

  # band
  d="$s/24_band_nosoc_HR"
  if [ -d "$d" ]; then
    if [ -s "$d/EIGENVAL" ] && [ -s "$d/OUTCAR" ]; then
      gap_line=$(python3 - "$d" <<'PY'
import sys, re, math
d=sys.argv[1]

# read E_F from OUTCAR
E_F=None
with open(f"{d}/OUTCAR") as f:
    for ln in f:
        if "E-fermi" in ln:
            m=re.search(r"E-fermi\s*:\s*([-\d\.Ee+]+)", ln)
            if m: E_F=float(m.group(1))

# read energies from EIGENVAL
E=[]
with open(f"{d}/EIGENVAL") as f:
    for ln in f:
        sp=ln.split()
        if len(sp)==3 and sp[0].lstrip("+-").isdigit():
            try: E.append(float(sp[1]))
            except: pass

# infer NBANDS from first k-block (optional, for direct/indirect)
nb=None
with open(f"{d}/EIGENVAL") as f:
    lines=f.readlines()
L=len(lines)
i=0
while i<L:
    sp=lines[i].split()
    if len(sp)==4:  # kx ky kz weight
        j=i+1; cnt=0
        while j<L:
            sp2=lines[j].split()
            if len(sp2)==3 and sp2[0].lstrip("+-").isdigit():
                cnt+=1; j+=1
            else:
                break
        if cnt>0: nb=cnt; break
    i+=1

# compute gap relative to E_F
vb = max([e for e in E if E_F is not None and e < E_F], default=math.nan)
cb = min([e for e in E if E_F is not None and e >= E_F], default=math.nan)
Eg = cb - vb if (not math.isnan(vb) and not math.isnan(cb)) else math.nan

# direct/indirect?
gtype="unknown"
if nb and not math.isnan(Eg):
    try:
        iv=E.index(vb); ic=E.index(cb)
        gtype="direct" if (iv//nb)==(ic//nb) else "indirect"
    except: pass

# NKPTS from OUTCAR (reliable)
NK="NA"
with open(f"{d}/OUTCAR") as f:
    for ln in f:
        if "NKPTS" in ln:
            m=re.search(r"NKPTS\s*=\s*(\d+)", ln)
            if m: NK=m.group(1); break

print(f"{Eg if not math.isnan(Eg) else float('nan'):.3f},{gtype},{NK}")
PY
)
      Eg=${gap_line%%,*}; rest=${gap_line#*,}
      gtype=${rest%%,*}; nk=${rest##*,}
      echo "$s, band, -, ${nk:-NA}, -, ${Eg:-NA}, ${gtype:-NA}"
    else
      echo "$s, band, -, NA, -, NA, NA"
    fi
  fi
done
