"""
Origin-style high-quality plots for publication/presentation
"""
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from pymatgen.io.vasp.outputs import Vasprun, Procar, Eigenval
from pymatgen.electronic_structure.core import Spin
from scipy.ndimage import gaussian_filter1d

# =====================================================
# Origin Style Configuration
# =====================================================
def setup_origin_style():
    """Configure matplotlib to mimic Origin software style"""

    # Font settings - Origin typically uses Arial
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
    rcParams['font.size'] = 14
    rcParams['axes.labelsize'] = 16
    rcParams['axes.titlesize'] = 16
    rcParams['xtick.labelsize'] = 14
    rcParams['ytick.labelsize'] = 14
    rcParams['legend.fontsize'] = 13

    # Line widths
    rcParams['axes.linewidth'] = 1.5
    rcParams['lines.linewidth'] = 2.0
    rcParams['xtick.major.width'] = 1.5
    rcParams['ytick.major.width'] = 1.5
    rcParams['xtick.minor.width'] = 1.0
    rcParams['ytick.minor.width'] = 1.0

    # Tick settings - Origin style with ticks inside
    rcParams['xtick.major.size'] = 6
    rcParams['ytick.major.size'] = 6
    rcParams['xtick.minor.size'] = 3
    rcParams['ytick.minor.size'] = 3
    rcParams['xtick.direction'] = 'in'
    rcParams['ytick.direction'] = 'in'
    rcParams['xtick.top'] = True
    rcParams['ytick.right'] = True

    # Grid (off by default in Origin)
    rcParams['axes.grid'] = False

    # Legend
    rcParams['legend.frameon'] = True
    rcParams['legend.fancybox'] = False
    rcParams['legend.edgecolor'] = 'black'
    rcParams['legend.shadow'] = False

    # Figure
    rcParams['figure.facecolor'] = 'white'
    rcParams['axes.facecolor'] = 'white'

    # Math text
    rcParams['mathtext.default'] = 'regular'

def get_origin_colors():
    """Return Origin-style color palette"""
    # Professional colors similar to Origin's default palette
    return [
        '#000000',  # Black
        '#FF0000',  # Red
        '#0000FF',  # Blue
        '#008000',  # Green
        '#FF00FF',  # Magenta
        '#FFA500',  # Orange
        '#800080',  # Purple
        '#00FFFF',  # Cyan
    ]

# =====================================================
# Data Reading Functions
# =====================================================
def read_doscar(doscar_file):
    """Read DOS from DOSCAR file"""
    with open(doscar_file, 'r') as f:
        lines = f.readlines()

    data = []
    for line in lines[6:]:
        parts = line.split()
        if len(parts) < 3:
            break
        E, DOS, iDOS = map(float, parts[:3])
        data.append([E, DOS])

    data = np.array(data)
    return data[:,0], data[:,1]

def process_strain_data(strain_dirs):
    """Process all strain directories and extract PDOS data"""
    results = {}

    for sd in strain_dirs:
        print(f"\nProcessing {sd}...")
        results[sd] = {}

        # Process SOC calculation for PDOS
        calc = "soc"
        vasprun_file = os.path.join(sd, calc, "vasprun.xml")
        eigen_file = os.path.join(sd, calc, "EIGENVAL")
        procar_file = os.path.join(sd, calc, "PROCAR")

        if not all(os.path.exists(f) for f in [vasprun_file, eigen_file, procar_file]):
            print(f"  [!] Missing required files in {sd}/{calc}")
            continue

        # Read vasprun for structure and Fermi energy
        vasprun = Vasprun(vasprun_file, parse_potcar_file=False)
        structure = vasprun.final_structure
        efermi = vasprun.efermi

        # Read eigenvalues
        eigen = Eigenval(eigen_file)
        eig_raw = eigen.eigenvalues

        try:
            occ_raw = eigen.occupancies
            if isinstance(eig_raw, dict):
                spin_key = list(eig_raw.keys())[0]
                energies = eig_raw[spin_key][:,:,0] - efermi
                occs = occ_raw[spin_key][:,:,0]
            else:
                energies = eig_raw[:,:,0] - efermi
                occs = occ_raw[:,:,0]
        except AttributeError:
            if isinstance(eig_raw, dict):
                spin_key = list(eig_raw.keys())[0]
                energies = eig_raw[spin_key][:,:,0] - efermi
                occs = eig_raw[spin_key][:,:,1]
            else:
                energies = eig_raw[:,:,0] - efermi
                occs = eig_raw[:,:,1]

        nk, nb = energies.shape

        # Calculate band edges
        occupied = energies[occs > 1e-3]
        empty = energies[occs < 1e-3]
        VBM = occupied.max()
        CBM = empty.min()
        gap = CBM - VBM

        print(f"  VBM={VBM:.3f} eV, CBM={CBM:.3f} eV, gap={gap:.3f} eV")

        # Read PROCAR for PDOS
        procar = Procar(procar_file)

        if isinstance(procar.data, dict):
            if Spin.up in procar.data:
                proj = procar.data[Spin.up]
            else:
                proj = list(procar.data.values())[0]
        else:
            proj = procar.data

        nk2, nb2, nions, norb = proj.shape
        nk_use = min(nk, nk2)
        nb_use = min(nb, nb2)

        # Find Mo atoms
        mo_idx = [i for i, s in enumerate(structure.sites) if s.species_string.lower() == "mo"]
        if len(mo_idx) == 0:
            mo_idx = [0]

        # Calculate PDOS for Mo d_z2 (orbital index 6)
        emin, emax = -3, 3
        bins = np.linspace(emin, emax, 600)
        dz2 = np.zeros_like(bins)

        for kpt in range(nk_use):
            for band in range(nb_use):
                e = energies[kpt, band]
                if emin < e < emax:
                    idx = np.searchsorted(bins, e) - 1
                    if 0 <= idx < len(bins):
                        dz2[idx] += proj[kpt, band, mo_idx, 6].sum()

        # Smooth the PDOS
        dz2 = gaussian_filter1d(dz2, sigma=2)

        results[sd] = {
            'bins': bins,
            'dz2': dz2,
            'VBM': VBM,
            'CBM': CBM,
            'gap': gap
        }

    return results

# =====================================================
# Plotting Functions
# =====================================================
def plot_pdos_comparison(results, strain_dirs, output_file="PDOS_comparison_origin.png"):
    """
    Create Origin-style PDOS comparison plot
    """
    setup_origin_style()
    colors = get_origin_colors()

    fig, ax = plt.subplots(figsize=(9, 6), dpi=600)

    # Strain labels for legend
    strain_labels = {
        'strain099': '-1% strain',
        'strain100': '0% strain',
        'strain101': '+1% strain'
    }

    # Plot each strain's PDOS
    for idx, sd in enumerate(strain_dirs):
        if sd in results and results[sd]['dz2'] is not None:
            label = strain_labels.get(sd, sd)
            ax.plot(results[sd]['bins'], results[sd]['dz2'],
                   color=colors[idx+1], linewidth=2.5, label=label)

    # Add Fermi level line
    ax.axvline(0, color='black', linestyle='--', linewidth=1.5,
              label='$E_F$', alpha=0.7)

    # Labels and title
    ax.set_xlabel('$E - E_F$ (eV)', fontsize=18, fontweight='bold')
    ax.set_ylabel('PDOS (Mo $d_{z^2}$)', fontsize=18, fontweight='bold')
    ax.set_xlim(-3, 3)

    # Set y-axis to start from 0
    ymin, ymax = ax.get_ylim()
    ax.set_ylim(0, ymax * 1.05)

    # Legend
    ax.legend(loc='upper right', frameon=True, edgecolor='black',
             fancybox=False, framealpha=1.0, fontsize=14)

    # Minor ticks
    ax.minorticks_on()

    # Tight layout
    plt.tight_layout()

    # Save with high DPI
    plt.savefig(output_file, dpi=600, bbox_inches='tight',
               facecolor='white', edgecolor='none')
    print(f"\n[+] Saved: {output_file}")

    plt.close()

def plot_pdos_with_bandedges(results, strain_dirs, output_file="PDOS_with_bands_origin.png"):
    """
    Create Origin-style PDOS comparison with VBM/CBM markers
    """
    setup_origin_style()
    colors = get_origin_colors()

    fig, ax = plt.subplots(figsize=(9, 6), dpi=600)

    strain_labels = {
        'strain099': '-1% strain',
        'strain100': '0% strain',
        'strain101': '+1% strain'
    }

    # Plot PDOS
    for idx, sd in enumerate(strain_dirs):
        if sd in results and results[sd]['dz2'] is not None:
            label = strain_labels.get(sd, sd)
            ax.plot(results[sd]['bins'], results[sd]['dz2'],
                   color=colors[idx+1], linewidth=2.5, label=label)

    # Add Fermi level
    ax.axvline(0, color='black', linestyle='--', linewidth=1.5,
              label='$E_F$', alpha=0.7)

    # Add VBM/CBM regions with shading
    ymin, ymax = ax.get_ylim()

    # Find average VBM and CBM
    avg_vbm = np.mean([results[sd]['VBM'] for sd in strain_dirs if sd in results])
    avg_cbm = np.mean([results[sd]['CBM'] for sd in strain_dirs if sd in results])

    # Shade band gap region
    ax.axvspan(avg_vbm, avg_cbm, alpha=0.15, color='gray', label='Band gap')

    # Labels
    ax.set_xlabel('$E - E_F$ (eV)', fontsize=18, fontweight='bold')
    ax.set_ylabel('PDOS (Mo $d_{z^2}$)', fontsize=18, fontweight='bold')
    ax.set_xlim(-3, 3)
    ax.set_ylim(0, ymax * 1.05)

    # Legend
    ax.legend(loc='upper right', frameon=True, edgecolor='black',
             fancybox=False, framealpha=1.0, fontsize=13)

    # Minor ticks
    ax.minorticks_on()

    plt.tight_layout()
    plt.savefig(output_file, dpi=600, bbox_inches='tight',
               facecolor='white', edgecolor='none')
    print(f"[+] Saved: {output_file}")

    plt.close()

# =====================================================
# Main Execution
# =====================================================
if __name__ == "__main__":
    print("="*60)
    print("Creating Origin-style high-quality plots")
    print("="*60)

    strain_dirs = ["strain099", "strain100", "strain101"]

    # Process data
    print("\nReading and processing data...")
    results = process_strain_data(strain_dirs)

    # Create plots
    print("\nGenerating plots...")
    plot_pdos_comparison(results, strain_dirs, "PDOS_comparison_origin.png")
    plot_pdos_with_bandedges(results, strain_dirs, "PDOS_with_bands_origin.png")

    print("\n" + "="*60)
    print("All Origin-style plots completed!")
    print("="*60)
    print("\nGenerated files:")
    print("  - PDOS_comparison_origin.png (600 dpi)")
    print("  - PDOS_with_bands_origin.png (600 dpi)")
    print("\nThese high-resolution images are ready for publication/presentation.")
