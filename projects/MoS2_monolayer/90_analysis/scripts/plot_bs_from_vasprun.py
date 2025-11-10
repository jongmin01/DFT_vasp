# Plot band structure from a *band* run vasprun.xml (line-mode KPOINTS).
# BSPlotter.get_plot() does NOT take an 'ax' kwarg; pass ylim instead.

from pathlib import Path
from pymatgen.io.vasp.outputs import BSVasprun
from pymatgen.electronic_structure.plotter import BSPlotter

# >>> Change this to your band directory (must contain vasprun.xml from BAND run)
band_dir = Path("strain_100/24_band_soc_HR")
vxml = band_dir / "vasprun.xml"

if not vxml.exists():
    raise FileNotFoundError(f"Band vasprun.xml not found at: {vxml}")

# Parse as symmetry-line band structure
bsv = BSVasprun(str(vxml), parse_projected_eigen=False)
bs = bsv.get_band_structure(line_mode=True)  # ensures BandStructureSymmLine

# Get a matplotlib Figure directly; cannot pass 'ax'
fig = BSPlotter(bs).get_plot(ylim=(-3, 3))

# Save and/or show
fig.savefig(band_dir / "band_soc.png", dpi=300, bbox_inches="tight")
print(f"[OK] Saved: {band_dir/'band_soc.png'}")
