#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plot_band_from_vasprun_v2_6.py

- Input: multiple vasprun.xml paths
- Output: PNG band plots saved in plots/bands/
- sumo NOT required (pymatgen + matplotlib only)
- Default behavior (no extra options needed):
    * orbital-projected colors (if available)
    * label edges (gap/type/vbm/cbm)
    * auto-detect SOC/NoSOC tag
    * auto-detect strain number from folder name

If projected data is missing → automatically falls back to normal band colors.
"""

import os
import re
import argparse
import matplotlib.pyplot as plt
from pymatgen.io.vasp.outputs import BSVasprun
from pymatgen.electronic_structure.bandstructure import BandStructureSymmLine
from pymatgen.electronic_structure.plotter import BSPlotter, BSPlotterProjected


def extract_strain(path):
    """Return strain as decimal string, e.g. strain_099 → 0.99."""
    m = re.search(r"strain[_-]?(\d+)", path.lower())
    if not m:
        return ""
    raw = m.group(1)
    try:
        if len(raw) == 3:
            return f"{int(raw)/100:.2f}"
        elif len(raw) == 2:
            return f"{int(raw)/10:.1f}"
        else:
            return raw
    except:
        return raw


def detect_tag(path):
    """Return SOC or NoSOC if detected in path, else ''."""
    p = path.lower()
    if "soc" in p and "nosoc" not in p:
        return "SOC"
    elif "nosoc" in p:
        return "NoSOC"
    return ""


def plot_single_band(xml_path, out_path,
                     emin, emax,
                     use_project, label_edges,
                     auto_tag, auto_strain):
    """Load vasprun.xml → shift VBM to 0 → plot → PNG."""
    print(f"[DEBUG] Reading: {xml_path}")

    # Load band structure
    v = BSVasprun(xml_path, parse_projected_eigen=use_project)
    bs = v.get_band_structure(line_mode=True)

    gap   = bs.get_band_gap()
    vbm   = bs.get_vbm()
    cbm   = bs.get_cbm()
    Eg    = float(gap.get("energy", 0.0))
    gtype = "direct" if gap.get("direct", False) else "indirect"

    vbm_e = float(vbm["energy"])

    # Force VBM = 0 eV by rewriting efermi
    d = bs.as_dict()
    d["efermi"] = vbm_e
    bs_vbm0 = BandStructureSymmLine.from_dict(d)

    # Auto tag
    tag = detect_tag(xml_path) if auto_tag else ""

    # Auto strain
    strain_txt = extract_strain(xml_path) if auto_strain else ""

    # Title
    title_txt = "Band Structure"
    if strain_txt:
        title_txt += f" (strain {strain_txt})"
    if tag:
        title_txt += f" – {tag}"

    # Try projected plot, fallback if no projected data
    if use_project:
        try:
            print("   -> Attempting orbital-projected plotting...")
            bsplt = BSPlotterProjected(bs_vbm0)
            projected_on = True
        except ValueError:
            print("   -> No projected data found. Falling back to normal band colors.")
            bsplt = BSPlotter(bs_vbm0)
            projected_on = False
    else:
        bsplt = BSPlotter(bs_vbm0)
        projected_on = False

    # Draw plot
    ax = bsplt.get_plot()
    ax.set_ylim(emin, emax)
    ax.axhline(0.0, ls="--", lw=0.8, alpha=0.6, color="black", label="VBM=0 eV")
    ax.set_ylabel("Energy (eV)")
    ax.set_title(title_txt)

    # Gap annotation
    if label_edges:
        def klabel(edge):
            kp = edge.get("kpoint", None)
            if kp is not None and getattr(kp, "label", None):
                return kp.label
            idx = edge.get("kpoint_index", None)
            return f"k#{idx}" if idx is not None else "unknown"

        info = f"Gap={Eg:.3f} eV ({gtype}), VBM@{klabel(vbm)}, CBM@{klabel(cbm)}"
        ax.text(0.02, 0.98, info, transform=ax.transAxes,
                va="top", ha="left", fontsize=9,
                bbox=dict(facecolor="white", alpha=0.6, edgecolor="none"))

    # Legend
    ax.legend(loc="lower left", fontsize=8)

    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    plt.close()
    print(f"[OK] Saved: {out_path}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-f", "--files", nargs="+", required=True,
                    help="Paths to vasprun.xml files")

    ap.add_argument("--emin", type=float, default=-3.0)
    ap.add_argument("--emax", type=float, default= 3.0)

    # Default settings: ON
    ap.add_argument("--no-project", action="store_true")
    ap.add_argument("--no-label",   action="store_true")
    ap.add_argument("--no-tag",     action="store_true")
    ap.add_argument("--no-strain",  action="store_true")

    args = ap.parse_args()

    # Defaults (unless disabled)
    use_project = not args.no_project
    label_edges = not args.no_label
    auto_tag    = not args.no_tag
    auto_strain = not args.no_strain

    out_dir = "plots/bands"
    os.makedirs(out_dir, exist_ok=True)

    for xml in args.files:
        if not os.path.isfile(xml):
            print(f"[WARN] File not found: {xml}")
            continue

        base = xml.replace("/", "_").replace("vasprun.xml", "")
        base = base.strip("_")
        out_path = os.path.join(out_dir, f"{base}.png")

        plot_single_band(
            xml_path=xml,
            out_path=out_path,
            emin=args.emin,
            emax=args.emax,
            use_project=use_project,
            label_edges=label_edges,
            auto_tag=auto_tag,
            auto_strain=auto_strain
        )

    print("\n[DONE] All band plots saved in plots/bands/")


if __name__ == "__main__":
    main()

