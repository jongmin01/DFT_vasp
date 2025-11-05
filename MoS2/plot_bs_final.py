#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Plot band structure + spin-resolved DOS (if SOC) side-by-side.
Works even if the bandstructure has no ticks (non-line-mode).
Output saved to plot/bands/.
"""

import os
import argparse
import numpy as np
import matplotlib.pyplot as plt

from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.electronic_structure.core import Spin
from pymatgen.electronic_structure.bandstructure import BandStructure

def short_name(xml_path):
    import re
    lower = xml_path.lower()
    tag = ""

    m = re.search(r"strain[_-]?0*([0-9]+)", lower)
    if m:
        tag += m.group(1)

    if "soc" in lower:
        tag += "_soc"
    elif ("nosoc" in lower) or ("no_soc" in lower):
        tag += "_nosoc"

    if tag == "":
        tag = "plot"
    return tag


def plot_bs_dos(xml_file, emin=-3, emax=3):

    print(f"\n▶ Processing : {xml_file}")

    if not os.path.exists(xml_file):
        print(f"❌ Cannot find: {xml_file}")
        return

    vr = Vasprun(xml_file, parse_projected_eigen=False, parse_potcar_file=False)
    bs = vr.get_band_structure()

    # Must be a band structure object
    if not isinstance(bs, BandStructure):
        print("❌ This is not a valid band structure")
        return

    # energy shift
    vbm_info = bs.get_vbm()
    cbm_info = bs.get_cbm()
    vbm_e = vbm_info["energy"]
    cbm_e = cbm_info["energy"]
    Eg = abs(cbm_e - vbm_e)

    bands = bs.bands

    # distances along k-path
    # if no ticks, distances still exist
    dist = np.array(bs.distance)

    # ---- Begin Plot ----
    fig = plt.figure(figsize=(9, 5))
    ax1 = fig.add_axes([0.1, 0.15, 0.55, 0.8])

    # plot each spin & band
    for spin in bands:
        eig = np.array(bands[spin]) - vbm_e
        for b in range(eig.shape[0]):
            ax1.plot(dist, eig[b], lw=1.3)

    ax1.set_ylim(emin, emax)
    ax1.axhline(0, color="gray", ls="--")
    ax1.set_ylabel("Energy (eV)")

    # Try ticks only if available
    try:
        ticks = bs.as_dict()["ticks"]
        xpos = ticks["distance"]
        labels = [x["label"] for x in bs.as_dict()["ticks"]["label"]]
        for t in xpos:
            ax1.axvline(t, color="black", lw=1)
        ax1.set_xticks(xpos)
        ax1.set_xticklabels(labels)
    except KeyError:
        # no ticks → no labels
        ax1.set_xticks([])

    ax1.set_title(f"Band structure (Eg={Eg:.3f} eV)")

    # mark VBM/CBM
    v_idx = vbm_info["kpoint_index"][0]
    c_idx = cbm_info["kpoint_index"][0]
    ax1.scatter(dist[v_idx], 0, color="red", zorder=5, label="VBM")
    ax1.scatter(dist[c_idx], cbm_e - vbm_e, color="green", zorder=5, label="CBM")
    ax1.legend(fontsize=8)

    # ---- DOS Panel ----
    ax2 = fig.add_axes([0.7, 0.15, 0.25, 0.8])

    doscar = os.path.join(os.path.dirname(xml_file), "DOSCAR")
    if os.path.exists(doscar):
        dos = vr.complete_dos
        energies = dos.energies - vbm_e
        try:
            up = dos.get_densities(spin=Spin.up)
            down = dos.get_densities(spin=Spin.down)
            ax2.plot(up, energies, color="blue", label="↑")
            ax2.plot(down, energies, color="orange", label="↓")
            ax2.legend(fontsize=8)
        except Exception:
            td = dos.get_densities()
            ax2.plot(td, energies, color="black")

        ax2.set_ylim(emin, emax)
        ax2.axhline(0, color="gray", ls="--")
        ax2.set_xlabel("DOS")
    else:
        ax2.text(0.2, 0.5, "No DOS", transform=ax2.transAxes)
        ax2.set_xticks([])
        ax2.set_yticks([])

    # Save
    out_dir = "plot/bands"
    os.makedirs(out_dir, exist_ok=True)
    name = short_name(xml_file)
    out = os.path.join(out_dir, f"{name}.png")
    plt.savefig(out, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"✅ Saved: {out}")


def main():
    parser = argparse.ArgumentParser(description="Plot BS + DOS")
    parser.add_argument("xml", nargs="+")
    parser.add_argument("--emin", type=float, default=-3)
    parser.add_argument("--emax", type=float, default=3)
    args = parser.parse_args()

    for f in args.xml:
        plot_bs_dos(f, args.emin, args.emax)


if __name__ == "__main__":
    main()

