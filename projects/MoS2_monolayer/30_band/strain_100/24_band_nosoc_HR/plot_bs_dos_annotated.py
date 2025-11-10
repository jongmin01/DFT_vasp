# file: plot_bs_dos_annotated.py
# Goal: Plot band (line-mode, SOC) + total DOS and annotate VBM/CBM, Eg, and high-symmetry ticks.
# Works across pymatgen versions by avoiding projection-only helpers.

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

from pymatgen.io.vasp.outputs import BSVasprun, Vasprun
from pymatgen.electronic_structure.core import Spin
from pymatgen.electronic_structure.plotter import BSPlotter

# ==== USER CONFIG ====
band_dir = Path("")   # <- BAND run (line-mode KPOINTS) with vasprun.xml
dos_dir  = band_dir                            # <- DOS-capable vasprun.xml (can be same dir)
ylim = (-3, 3)                                 # y-range around EF for both panels

# ==== Load inputs ====
vxml_band = band_dir / "vasprun.xml"
vxml_dos  = dos_dir  / "vasprun.xml"
if not vxml_band.exists():
    raise FileNotFoundError(f"Band vasprun.xml not found: {vxml_band}")
if not vxml_dos.exists():
    raise FileNotFoundError(f"DOS vasprun.xml not found:  {vxml_dos}")

# Parse band (symmetry lines)
bsv = BSVasprun(str(vxml_band), parse_projected_eigen=False)
bs  = bsv.get_band_structure(line_mode=True)  # BandStructureSymmLine

# Extract distances (x), ticks (vertical lines) and labels
def get_ticks_safe(bsobj):
    """Return (tick_positions, tick_labels). Falls back to k-point labels if needed."""
    try:
        ticks = bsobj.get_ticks()
        return ticks["distance"], ticks["label"]
    except Exception:
        # Fallback: collect k-points with a label and use their distances
        dists, labels = [], []
        for kp, dist in zip(bsobj.kpoints, bsobj.distance):
            if kp.label:
                dists.append(dist)
                labels.append(kp.label)
        # Ensure at least segment boundaries (Γ, K, M, Γ) appear
        return dists, labels

tick_pos, tick_lbl = get_ticks_safe(bs)

# Energies array and distances
efermi = bs.efermi
# Choose first available spin
spin = list(bs.bands.keys())[0]
E = np.array(bs.bands[spin])  # shape: (nbands, nkpts)
# Shift energies to EF=0 for plotting and gap search
E_shift = E - efermi
x = np.array(bs.distance)

# Compute VBM/CBM and gap
def vbm_cbm_from_bands(E_shift):
    """Return (VBM_e, VBM_k, VBM_band), (CBM_e, CBM_k, CBM_band), Eg"""
    nb, nk = E_shift.shape
    vbm_e, vbm_k, vbm_b = -1e9, -1, -1
    cbm_e, cbm_k, cbm_b = +1e9, -1, -1
    for k in range(nk):
        col = E_shift[:, k]
        # occupied bands: e <= 0; empty: e > 0
        occ_idx = np.where(col <= 0.0)[0]
        emp_idx = np.where(col >  0.0)[0]
        if occ_idx.size > 0:
            e_top = col[occ_idx].max()
            b_top = occ_idx[col[occ_idx].argmax()]
            if e_top > vbm_e:
                vbm_e, vbm_k, vbm_b = float(e_top), k, int(b_top)
        if emp_idx.size > 0:
            e_bot = col[emp_idx].min()
            b_bot = emp_idx[col[emp_idx].argmin()]
            if e_bot < cbm_e:
                cbm_e, cbm_k, cbm_b = float(e_bot), k, int(b_bot)
    Eg = max(0.0, cbm_e - vbm_e)
    return (vbm_e, vbm_k, vbm_b), (cbm_e, cbm_k, cbm_b), Eg

(v_e, v_k, v_b), (c_e, c_k, c_b), Eg = vbm_cbm_from_bands(E_shift)
is_direct = (v_k == c_k)

# Prepare DOS
vr_dos = Vasprun(str(vxml_dos), parse_projected_eigen=False)
cdos   = vr_dos.complete_dos
common_efermi = bs.efermi
E_shift = E - common_efermi
E_dos   = cdos.energies - common_efermi
D_up   = cdos.densities.get(Spin.up, None)
D_dn   = cdos.densities.get(Spin.down, None)

# ==== Plot ====
fig, (ax_b, ax_d) = plt.subplots(1, 2, figsize=(10, 4), gridspec_kw={"width_ratios": [2.2, 1.0]})

# Left: band structure
for b in range(E_shift.shape[0]):
    ax_b.plot(x, E_shift[b, :], lw=1.0)

# High-symmetry vertical lines and labels
for t in tick_pos:
    ax_b.axvline(t, color="k", lw=1.8)
# If matplotlib placed labels from BSPlotter, great. Otherwise, add simple ticks for Γ–K–M–Γ
# Try to use BSPlotter just to draw its x-ticks nicely, then copy them.
try:
    bsplt = BSPlotter(bs)
    tmp = bsplt.get_plot(ylim=ylim)
    # figure or axes, we only need xticks/labels
    ax_src = tmp.axes[0] if hasattr(tmp, "axes") else tmp
    ax_b.set_xticks(ax_src.get_xticks())
    ax_b.set_xticklabels([lbl.get_text() for lbl in ax_src.get_xticklabels()])
    plt.close(tmp if hasattr(tmp, "savefig") else None)
except Exception:
    if tick_pos and tick_lbl:
        ax_b.set_xticks(tick_pos)
        ax_b.set_xticklabels(tick_lbl)

ax_b.axhline(0.0, ls="--", lw=0.8, color="tab:blue", alpha=0.8)
ax_b.set_title("Band structure (SOC)")
ax_b.set_ylim(*ylim)
ax_b.set_ylabel("Energy (eV)")

# Annotate VBM/CBM and Eg
ax_b.scatter([x[v_k]], [v_e], color="red", zorder=5, label="VBM")
ax_b.scatter([x[c_k]], [c_e], color="green", zorder=5, label="CBM")
txt = f"Eg = {Eg:.3f} eV  ({'direct' if is_direct else 'indirect'})"
# Place text near the top-left
ax_b.text(0.02, 0.95, txt, transform=ax_b.transAxes, va="top", ha="left",
          bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="gray", alpha=0.8), fontsize=10)
ax_b.legend(frameon=False, loc="lower left")

# Right: total DOS
if D_up is not None:
    ax_d.plot(D_up, E_dos, lw=1.0)
if D_dn is not None:
    ax_d.plot(-D_dn, E_dos, lw=1.0)  # mirror spin-down if present
ax_d.axhline(0.0, ls="--", lw=0.8, color="tab:blue", alpha=0.8)
ax_d.set_ylim(*ylim)
ax_d.set_xlabel("DOS")
ax_d.set_yticklabels([])  # keep energy only on band panel

plt.tight_layout()
out = band_dir / "band_dos_annotated.png"
fig.savefig(out, dpi=300, bbox_inches="tight")
print(f"[OK] Saved: {out}")
