# file: plot_band_and_total_dos_fixed.py
# Robust: works whether BSPlotter.get_plot() returns a Figure or an Axes.
# Plots band (from line-mode vasprun.xml) + total DOS (from any vasprun.xml).

from pathlib import Path
import matplotlib.pyplot as plt
from pymatgen.io.vasp.outputs import BSVasprun, Vasprun
from pymatgen.electronic_structure.plotter import BSPlotter
from pymatgen.electronic_structure.core import Spin

# ===== PATHS (edit) =====
band_dir = Path("strain_100/24_band_soc_HR")   # band run (line-mode KPOINTS)
dos_dir  = band_dir                            # or a separate DOS run dir

vxml_band = band_dir / "vasprun.xml"
vxml_dos  = dos_dir  / "vasprun.xml"
if not vxml_band.exists():
    raise FileNotFoundError(f"Band vasprun.xml not found: {vxml_band}")
if not vxml_dos.exists():
    raise FileNotFoundError(f"DOS vasprun.xml not found:  {vxml_dos}")

# ===== Load band (symmetry lines) =====
bsv = BSVasprun(str(vxml_band), parse_projected_eigen=False)
bs  = bsv.get_band_structure(line_mode=True)  # BandStructureSymmLine

# ===== Load total DOS =====
vr_dos = Vasprun(str(vxml_dos), parse_projected_eigen=False)
cdos   = vr_dos.complete_dos
energies = cdos.energies - cdos.efermi
dos_up   = cdos.densities.get(Spin.up, None)
dos_dn   = cdos.densities.get(Spin.down, None)

# ===== Figure =====
fig, (ax_band, ax_dos) = plt.subplots(1, 2, figsize=(9, 4), gridspec_kw={"width_ratios": [2, 1]})

# --- Band: draw using BSPlotter, then copy artists to our axes
bsplt = BSPlotter(bs)
maybe_fig = bsplt.get_plot(ylim=(-3, 3))  # returns Figure OR Axes depending on version

# Get the real Axes that contains the band lines
if isinstance(maybe_fig, plt.Figure):
    src_ax = maybe_fig.axes[0]
else:
    # already an Axes
    src_ax = maybe_fig

# copy line artists
for line in src_ax.get_lines():
    ax_band.plot(line.get_xdata(), line.get_ydata(),
                 color=line.get_color(), lw=line.get_linewidth())

# decorate band axes
ax_band.axhline(0.0, ls="--", lw=0.8)
ax_band.set_title("Band structure (SOC)")
ax_band.set_ylim(-3, 3)
ax_band.set_xlabel("")
ax_band.set_ylabel("Energy (eV)")

# Close temp Figure ONLY if it is a Figure
if isinstance(maybe_fig, plt.Figure):
    plt.close(maybe_fig)

# --- DOS: total DOS only (no projections needed)
if dos_up is not None:
    ax_dos.plot(dos_up, energies, label="↑")
if dos_dn is not None:
    ax_dos.plot(-dos_dn, energies, label="↓")
ax_dos.axhline(0.0, ls="--", lw=0.8)
ax_dos.set_ylim(-3, 3)
ax_dos.set_xlabel("DOS")
ax_dos.set_ylabel("")
if dos_dn is not None:
    ax_dos.legend(frameon=False)

plt.tight_layout()
out_png = band_dir / "band_and_total_dos_fixed.png"
fig.savefig(out_png, dpi=300, bbox_inches="tight")
print(f"[OK] Saved: {out_png}")
