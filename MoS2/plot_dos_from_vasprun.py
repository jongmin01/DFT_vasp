#!/usr/bin/env python3
import argparse, os
from pymatgen.io.vasp.outputs import Vasprun
import matplotlib.pyplot as plt
from pymatgen.electronic_structure.core import Spin

def _is_complete_xml(path: str) -> bool:
    try:
        with open(path, "rb") as f:
            f.seek(-256, os.SEEK_END)
            return b"</modeling>" in f.read()
    except Exception:
        return False

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-f", "--file", required=True)  # path to vasprun.xml
    ap.add_argument("-o", "--out", default="dos.png")
    ap.add_argument("--emin", type=float, default=-3.0)
    ap.add_argument("--emax", type=float, default= 3.0)
    ap.add_argument("--xmax", type=float, default=None)
    args = ap.parse_args()

    if not os.path.isfile(args.file) or not _is_complete_xml(args.file):
        raise SystemExit(f"[ERR] unusable vasprun.xml: {args.file}")

    v = Vasprun(args.file, parse_potcar_file=False)
    cdos = v.complete_dos

    e = cdos.energies - cdos.efermi
    plt.figure(figsize=(4,5), dpi=200)
    ax = plt.gca()

    for sp, dens in cdos.densities.items():
        x = dens if sp != Spin.down else -dens
        ax.plot(x, e, lw=1.3)

    ax.axhline(0.0, color="tab:blue", ls="--", lw=1)
    ax.set_ylim(args.emin, args.emax)
    if args.xmax:
        ax.set_xlim(-args.xmax, args.xmax)
    ax.set_xlabel("DOS")
    ax.set_ylabel("Energy (eV)")
    plt.tight_layout()
    plt.savefig(args.out)
    print(f"[OK] Saved {args.out}")

if __name__ == "__main__":
    main()
