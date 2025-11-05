#!/usr/bin/env python3
import os, re, glob, numpy as np

print("strain_dir,  SCF_converged,  TotalE(eV),  VBM(eV),  CBM(eV),  Gap(eV),  Type")

for sdir in sorted(glob.glob("strain_*")):
    scf = os.path.join(sdir, "12_scf_nosoc_hr", "OUTCAR")
    band = os.path.join(sdir, "24_band_nosoc_HR", "EIGENVAL")

    if not os.path.exists(scf) or not os.path.exists(band):
        print(f"{sdir:12s}  MISSING files")
        continue

    # --- SCF convergence ---
    converged = "reached" in open(scf).read()
    Etot = None
    with open(scf) as f:
        for line in f:
            if "free  energy   TOTEN" in line:
                Etot = float(line.split()[-2])
    conv_tag = "YES" if converged else "NO"

    # --- Band gap parsing from EIGENVAL ---
    with open(band) as f:
        lines = f.readlines()
    header = lines[7].split()
    nkpts, nbands = int(header[1]), int(header[2])
    E = []
    for ln in lines[8:]:
        parts = ln.split()
        if len(parts) == 0:
            continue
        if len(parts) == 4:  # k-point line
            continue
        if len(parts) == 3:  # band line
            E.append(float(parts[1]))
    E = np.array(E).reshape(nkpts, nbands)

    vbm = np.max(E[E < 0])
    cbm = np.min(E[E >= 0])
    gap = cbm - vbm
    Eg_type = "direct" if np.isclose(np.argmax(E==vbm)//nbands,
                                     np.argmax(E==cbm)//nbands, atol=1) else "indirect"

    print(f"{sdir:12s}  {conv_tag:6s}  {Etot:10.3f}  {vbm:8.3f}  {cbm:8.3f}  {gap:6.3f}  {Eg_type}")
