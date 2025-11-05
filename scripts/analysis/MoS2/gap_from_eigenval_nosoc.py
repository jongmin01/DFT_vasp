#!/usr/bin/env python3
# Robust band-gap from EIGENVAL (line-mode, NoSOC or SOC-off)
# - Detect k-block headers (4 floats)
# - Band lines: ISPIN=1 (3 cols) or ISPIN=2 (5 cols, unknown order)
# - Use occupation threshold (occ>0.5 occupied, occ<0.5 empty)
# - Report global VBM/CBM, Eg, direct/indirect, and k-index
# Optional: map k-index to Γ/K/M labels with --segments "80,80,80"

import sys, math, argparse

def is_float(s):
    try:
        float(s.replace('D','E'))
        return True
    except:
        return False

def split_e_occ_5cols(sp2):
    vals = [float(x) for x in sp2[1:5]]
    occ_idx = [i for i,v in enumerate(vals) if -1e-6 <= v <= 1.0+1e-6]
    eng_idx = [i for i in range(4) if i not in occ_idx]
    pairs = []
    if len(occ_idx)==2 and len(eng_idx)==2:
        # pair energy with nearest occupancy by index distance
        occs = occ_idx[:]
        for ei in eng_idx:
            oi = min(occs, key=lambda k: abs(k-ei))
            pairs.append((vals[ei], vals[oi]))
            occs.remove(oi)
    return pairs

def label_from_index(k, segments):
    if not segments: return f"k={k}"
    cum=[0]
    for s in segments: cum.append(cum[-1]+s)
    labs=["Γ","K","M","Γ"]
    for i in range(3):
        if cum[i] <= k < cum[i+1]:
            t = (k-cum[i])/(segments[i]-1) if segments[i]>1 else 0.0
            return f"{labs[i]}→{labs[i+1]} (idx {k}, t={t:.2f})"
    return f"k={k}"

def main():
    ap = argparse.ArgumentParser(description="Band gap from EIGENVAL (occupation-based, k-resolved).")
    ap.add_argument("-f","--file", default="EIGENVAL")
    ap.add_argument("--segments", help="comma-separated counts per segment (e.g. 80,80,80)", default=None)
    args = ap.parse_args()

    segments = None
    if args.segments:
        segments = [int(x) for x in args.segments.split(",")]

    with open(args.file) as f:
        lines=f.readlines()

    L=len(lines); i=0; kidx=-1
    VBM_k={}; CBM_k={}

    while i<L:
        sp=lines[i].split()
        if len(sp)==4 and all(is_float(x) for x in sp):  # k-header
            kidx+=1
            j=i+1
            occ_E=[]; unocc_E=[]
            while j<L:
                sp2=lines[j].split()
                if not sp2 or not sp2[0].lstrip("+-").isdigit():
                    break
                if len(sp2)==3 and all(is_float(x) for x in sp2[1:]):
                    E=float(sp2[1]); occ=float(sp2[2])
                    if occ>0.5: occ_E.append(E)
                    elif occ<0.5: unocc_E.append(E)
                    j+=1
                elif len(sp2)>=5 and all(is_float(x) for x in sp2[1:5]):
                    for E,occ in split_e_occ_5cols(sp2):
                        if 0.0<=occ<=1.0:
                            if occ>0.5: occ_E.append(E)
                            elif occ<0.5: unocc_E.append(E)
                    j+=1
                else:
                    break
            if occ_E: VBM_k[kidx]=max(occ_E)
            if unocc_E: CBM_k[kidx]=min(unocc_E)
            i=j
        else:
            i+=1

    if not VBM_k or not CBM_k:
        print(f"[ERR] Missing VBM/CBM at k-points (VBM_k={len(VBM_k)}, CBM_k={len(CBM_k)}). "
              "Increase NBANDS or check EIGENVAL format.")
        sys.exit(1)

    k_vbm, VBM = max(VBM_k.items(), key=lambda kv: kv[1])
    k_cbm, CBM = min(CBM_k.items(), key=lambda kv: kv[1])
    Eg = CBM - VBM
    gtype = "direct" if k_vbm==k_cbm else "indirect"

    lab_v = label_from_index(k_vbm, segments)
    lab_c = label_from_index(k_cbm, segments)

    print("=== Band gap (occupation-based, k-resolved) ===")
    print(f"num k with VBM_k : {len(VBM_k)}")
    print(f"num k with CBM_k : {len(CBM_k)}")
    print(f"VBM : {VBM:.6f} eV @ {lab_v}")
    print(f"CBM : {CBM:.6f} eV @ {lab_c}")
    print(f"Eg  : {Eg:.6f} eV ({gtype})")

if __name__ == "__main__":
    main()
