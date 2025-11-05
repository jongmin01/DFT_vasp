#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import argparse
import matplotlib.pyplot as plt
from pymatgen.io.vasp.outputs import BSVasprun
from pymatgen.electronic_structure.bandstructure import BandStructureSymmLine
from pymatgen.electronic_structure.plotter import BSPlotter

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-f", "--file", required=True, help="vasprun.xml (band run)")
    ap.add_argument("-o", "--out",  required=True, help="output image path (png/pdf)")
    ap.add_argument("--emin", type=float, default=-3.0)
    ap.add_argument("--emax", type=float, default= 3.0)
    ap.add_argument("--mark-edges", action="store_true")   # 예약(좌표 마킹은 라벨만)
    ap.add_argument("--label-edges", action="store_true")  # 텍스트 라벨 표시
    ap.add_argument("--marker-size", type=float, default=36.0)
    args = ap.parse_args()

    v = BSVasprun(args.file, parse_projected_eigen=False)
    bs = v.get_band_structure(line_mode=True)

    gap = bs.get_band_gap()
    vbm = bs.get_vbm()
    cbm = bs.get_cbm()
    vbm_e = vbm["energy"]
    cbm_e = cbm["energy"]
    Eg    = float(gap.get("energy", 0.0))
    gtype = "direct" if gap.get("direct", False) else "indirect"

    # VBM을 0 eV로 보이게: efermi를 VBM 에너지로 덮어씀
    d = bs.as_dict()
    d["efermi"] = vbm_e
    bs_vbm0: BandStructureSymmLine = BandStructureSymmLine.from_dict(d)

    bsplt = BSPlotter(bs_vbm0)
    ax = bsplt.get_plot()   # 구버전: Axes 반환
    ax.set_ylim(args.emin, args.emax)
    ax.axhline(0.0, ls="--", lw=0.8, alpha=0.6)  # VBM=0 기준선

    # 텍스트 라벨(안전)
    if args.label_edges:
        def klabel(edge):
            kp = edge.get("kpoint", None)
            if kp is not None and getattr(kp, "label", None):
                return kp.label
            idx = edge.get("kpoint_index", None)
            return f"k#{idx}" if idx is not None else "unknown"
        txt = f"Gap={Eg:.3f} eV ({gtype}) | VBM@{klabel(vbm)} | CBM@{klabel(cbm)}"
        ax.text(0.02, 0.98, txt, transform=ax.transAxes, va="top", ha="left", fontsize=9)

    plt.tight_layout()
    plt.savefig(args.out, dpi=300)
    plt.close()

if __name__ == "__main__":
    main()
