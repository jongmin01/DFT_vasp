"""
Comprehensive Origin-style high-quality plots for all analysis
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

def get_origin_markers():
    """Return Origin-style marker styles"""
    return ['o', 's', '^', 'D', 'v', '<', '>', 'p']

# =====================================================
# Data Reading Functions
# =====================================================
def read_doscar(doscar_file):
    """Read DOS from DOSCAR file"""
    with open(doscar_file, 'r') as f:
        lines = f.readlines()


    # Find where data actually starts
    for i, line in enumerate(lines):
        parts = line.split()
        if len(parts) >= 3:
            try:
                float(parts[0])
                start_idx = i
                break
            except ValueError:
                continue


    data = []
    for line in lines[6:]:
        parts = line.split()
        if len(parts) < 3:
            break
        E, DOS, iDOS = map(float, parts[:3])
        data.append([E, DOS])

    data = np.array(data)
    return data[:,0], data[:,1]

def process_all_data(strain_dirs, calc_types):
    """Process all strain directories and both SOC/noSOC calculations"""
    results = {}

    for sd in strain_dirs:
        print(f"\n{'='*50}")
        print(f"Processing {sd}")
        print('='*50)
        results[sd] = {}

        for calc in calc_types:
            print(f"\n  [{calc.upper()}]")
            vasprun_file = os.path.join(sd, calc, "vasprun.xml")
            eigen_file = os.path.join(sd, calc, "EIGENVAL")
            doscar_file = os.path.join(sd, calc, "DOSCAR")

            if not os.path.exists(vasprun_file) or not os.path.exists(eigen_file):
                print(f"    [!] Missing required files")
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

            print(f"    VBM={VBM:.3f} eV, CBM={CBM:.3f} eV, gap={gap:.3f} eV")

            # Read DOS
            DOS_bins, DOS_val = None, None
            if os.path.exists(doscar_file):
                Eb, Db = read_doscar(doscar_file)
                Eb -= efermi

                # Remove negative and very small values before smoothing
                Db[Db < 0] = 0

                # Smooth with smaller sigma
                DOS_val = gaussian_filter1d(Db, sigma=1)
                
                # ✅ Find near-EF region (±0.3 eV) and set baseline
                ef_window = (Eb > -0.3) & (Eb < 0.3)
                baseline = np.mean(DOS_val[ef_window])
                DOS_val -= baseline

                # ✅ Force any small negatives to 0
                DOS_val[DOS_val < 0] = 0

                DOS_bins = Eb

            # Initialize result dictionary for this calculation
            results[sd][calc] = {
                'VBM': VBM,
                'CBM': CBM,
                'gap': gap,
                'DOS_bins': DOS_bins,
                'DOS': DOS_val,
                'bins': None,
                'dz2': None
            }

            # For SOC, also calculate PDOS
            if calc == "soc":
                procar_file = os.path.join(sd, calc, "PROCAR")

                if os.path.exists(procar_file):
                    print(f"    Reading PROCAR for PDOS...")
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
                    mo_idx = [i for i, s in enumerate(structure.sites)
                             if s.species_string.lower() == "mo"]
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

                    results[sd][calc]['bins'] = bins
                    results[sd][calc]['dz2'] = dz2

    return results

# =====================================================
# Plotting Functions
# =====================================================

def plot_pdos_comparison(results, strain_dirs):
    """PDOS comparison across strains (SOC only)"""
    setup_origin_style()
    colors = get_origin_colors()

    fig, ax = plt.subplots(figsize=(9, 6), dpi=600)

    strain_labels = {
        'strain099': '-1% strain',
        'strain100': '0% strain',
        'strain101': '+1% strain'
    }

    for idx, sd in enumerate(strain_dirs):
        if sd in results and 'soc' in results[sd] and results[sd]['soc']['dz2'] is not None:
            label = strain_labels.get(sd, sd)
            ax.plot(results[sd]['soc']['bins'], results[sd]['soc']['dz2'],
                   color=colors[idx+1], linewidth=2.5, label=label)

    ax.axvline(0, color='black', linestyle='--', linewidth=1.5,
              label='$E_F$', alpha=0.7)

    ax.set_xlabel('$E - E_F$ (eV)', fontsize=18, fontweight='bold')
    ax.set_ylabel('PDOS (Mo $d_{z^2}$)', fontsize=18, fontweight='bold')
    ax.set_xlim(-3, 3)

    ymin, ymax = ax.get_ylim()
    ax.set_ylim(0, ymax * 1.05)

    ax.legend(loc='upper right', frameon=True, edgecolor='black',
             fancybox=False, framealpha=1.0, fontsize=14)

    ax.minorticks_on()
    plt.tight_layout()

    plt.savefig("origin_PDOS_comparison.png", dpi=600, bbox_inches='tight',
               facecolor='white', edgecolor='none')
    print("[+] Saved: origin_PDOS_comparison.png")
    plt.close()


def plot_dos_comparison_nosoc(results, strain_dirs):
    """DOS comparison for noSOC calculations"""
    setup_origin_style()
    colors = get_origin_colors()

    fig, ax = plt.subplots(figsize=(9, 6), dpi=600)

    strain_labels = {
        'strain099': '-1% strain',
        'strain100': '0% strain',
        'strain101': '+1% strain'
    }

    for idx, sd in enumerate(strain_dirs):
        if sd in results and 'nosoc' in results[sd] and results[sd]['nosoc']['DOS'] is not None:
            label = strain_labels.get(sd, sd)
            ax.plot(results[sd]['nosoc']['DOS_bins'], results[sd]['nosoc']['DOS'],
                   color=colors[idx+1], linewidth=2.5, label=label)

    ax.axvline(0, color='black', linestyle='--', linewidth=1.5,
              label='$E_F$', alpha=0.7)

    ax.set_xlabel('$E - E_F$ (eV)', fontsize=18, fontweight='bold')
    ax.set_ylabel('DOS (states/eV)', fontsize=18, fontweight='bold')
    ax.set_xlim(-3, 3)

    ymin, ymax = ax.get_ylim()
    #ax.set_ylim(0, ymax * 1.05)
    ax.set_ylim(-0.05, ymax * 1.05)

    ax.legend(loc='upper right', frameon=True, edgecolor='black',
             fancybox=False, framealpha=1.0, fontsize=14)

    ax.minorticks_on()
    plt.tight_layout()

    # Draw y=0 line in gray instead of black
    ax.axhline(0, color='gray', linewidth=1.0, alpha=0.4, zorder=1)

    # Bring plotted DOS lines to the very front
    for line in ax.lines:
        line.set_zorder(10)

    # Remove bottom spine overlap (optional)
    ax.spines['bottom'].set_color('gray')
    ax.spines['bottom'].set_linewidth(1.0)

    plt.savefig("origin_DOS_comparison_noSOC.png", dpi=600, bbox_inches='tight',
               facecolor='white', edgecolor='none')
    print("[+] Saved: origin_DOS_comparison_noSOC.png")
    plt.close()


def plot_soc_vs_nosoc_dos(results, strain_dirs):
    """SOC vs noSOC DOS comparison for each strain"""
    setup_origin_style()

    strain_labels = {
        'strain099': '-1% strain',
        'strain100': '0% strain',
        'strain101': '+1% strain'
    }

    for sd in strain_dirs:
        if sd not in results:
            continue

        if 'soc' not in results[sd] or 'nosoc' not in results[sd]:
            continue

        if results[sd]['soc']['DOS'] is None or results[sd]['nosoc']['DOS'] is None:
            continue

        fig, ax = plt.subplots(figsize=(9, 6), dpi=600)

        # Plot SOC
        ax.plot(results[sd]['soc']['DOS_bins'], results[sd]['soc']['DOS'],
               color='#FF0000', linewidth=2.5, label='SOC')

        # Plot noSOC
        ax.plot(results[sd]['nosoc']['DOS_bins'], results[sd]['nosoc']['DOS'],
               color='#0000FF', linewidth=2.5, linestyle='--', label='noSOC')

        ax.axvline(0, color='black', linestyle='--', linewidth=1.5,
                  label='$E_F$', alpha=0.7)

        ax.set_xlabel('$E - E_F$ (eV)', fontsize=18, fontweight='bold')
        ax.set_ylabel('DOS (states/eV)', fontsize=18, fontweight='bold')
        ax.set_xlim(-3, 3)

        ymin, ymax = ax.get_ylim()
        ax.set_ylim(0, ymax * 1.05)

        title = strain_labels.get(sd, sd)
        ax.set_title(title, fontsize=18, fontweight='bold', pad=15)

        ax.legend(loc='upper right', frameon=True, edgecolor='black',
                 fancybox=False, framealpha=1.0, fontsize=14)

        ax.minorticks_on()
        plt.tight_layout()

        filename = f"origin_DOS_SOC_vs_noSOC_{sd}.png"
        plt.savefig(filename, dpi=600, bbox_inches='tight',
                   facecolor='white', edgecolor='none')
        print(f"[+] Saved: {filename}")
        plt.close()


def plot_vbm_cbm_vs_strain(results, strain_dirs):
    """VBM and CBM vs strain plot"""
    setup_origin_style()
    colors = get_origin_colors()
    markers = get_origin_markers()

    fig, ax = plt.subplots(figsize=(9, 6), dpi=600)

    # Extract strain values
    strain_vals = []
    VBM_soc, CBM_soc = [], []
    VBM_nosoc, CBM_nosoc = [], []

    for sd in strain_dirs:
        sval = int(sd[-3:]) - 100
        strain_vals.append(sval)

        if 'soc' in results[sd]:
            VBM_soc.append(results[sd]['soc']['VBM'])
            CBM_soc.append(results[sd]['soc']['CBM'])
        else:
            VBM_soc.append(None)
            CBM_soc.append(None)

        if 'nosoc' in results[sd]:
            VBM_nosoc.append(results[sd]['nosoc']['VBM'])
            CBM_nosoc.append(results[sd]['nosoc']['CBM'])
        else:
            VBM_nosoc.append(None)
            CBM_nosoc.append(None)

    # Plot VBM
    ax.plot(strain_vals, VBM_soc, marker=markers[0], markersize=10,
           color=colors[1], linewidth=2.5, label='VBM (SOC)')
    ax.plot(strain_vals, VBM_nosoc, marker=markers[1], markersize=10,
           color=colors[1], linewidth=2.5, linestyle='--', label='VBM (noSOC)')

    # Plot CBM
    ax.plot(strain_vals, CBM_soc, marker=markers[0], markersize=10,
           color=colors[2], linewidth=2.5, label='CBM (SOC)')
    ax.plot(strain_vals, CBM_nosoc, marker=markers[1], markersize=10,
           color=colors[2], linewidth=2.5, linestyle='--', label='CBM (noSOC)')

    ax.set_xlabel('Strain (%)', fontsize=18, fontweight='bold')
    ax.set_ylabel('Energy (eV)', fontsize=18, fontweight='bold')

    ax.legend(loc='best', frameon=True, edgecolor='black',
             fancybox=False, framealpha=1.0, fontsize=13, ncol=2)

    ax.minorticks_on()
    plt.tight_layout()

    plt.savefig("origin_VBM_CBM_vs_strain.png", dpi=600, bbox_inches='tight',
               facecolor='white', edgecolor='none')
    print("[+] Saved: origin_VBM_CBM_vs_strain.png")
    plt.close()


def plot_bandgap_vs_strain(results, strain_dirs):
    """Bandgap vs strain plot"""
    setup_origin_style()
    colors = get_origin_colors()
    markers = get_origin_markers()

    fig, ax = plt.subplots(figsize=(9, 6), dpi=600)

    # Extract strain values
    strain_vals = []
    gap_soc, gap_nosoc = [], []

    for sd in strain_dirs:
        sval = int(sd[-3:]) - 100
        strain_vals.append(sval)

        if 'soc' in results[sd]:
            gap_soc.append(results[sd]['soc']['gap'])
        else:
            gap_soc.append(None)

        if 'nosoc' in results[sd]:
            gap_nosoc.append(results[sd]['nosoc']['gap'])
        else:
            gap_nosoc.append(None)

    # Plot bandgaps
    ax.plot(strain_vals, gap_soc, marker=markers[0], markersize=12,
           color=colors[2], linewidth=2.5, label='SOC')
    ax.plot(strain_vals, gap_nosoc, marker=markers[1], markersize=12,
           color=colors[1], linewidth=2.5, linestyle='--', label='noSOC')

    ax.set_xlabel('Strain (%)', fontsize=18, fontweight='bold')
    ax.set_ylabel('Band gap (eV)', fontsize=18, fontweight='bold')

    # Set y-axis to start near minimum value
    all_gaps = [g for g in gap_soc + gap_nosoc if g is not None]
    ymin = min(all_gaps) * 0.9
    ymax = max(all_gaps) * 1.1
    ax.set_ylim(ymin, ymax)

    ax.legend(loc='best', frameon=True, edgecolor='black',
             fancybox=False, framealpha=1.0, fontsize=14)

    ax.minorticks_on()
    plt.tight_layout()

    plt.savefig("origin_bandgap_vs_strain.png", dpi=600, bbox_inches='tight',
               facecolor='white', edgecolor='none')
    print("[+] Saved: origin_bandgap_vs_strain.png")
    plt.close()


# =====================================================
# Main Execution
# =====================================================
if __name__ == "__main__":
    print("\n" + "="*60)
    print(" Creating ALL Origin-style high-quality plots")
    print("="*60 + "\n")

    strain_dirs = ["strain099", "strain100", "strain101"]
    calc_types = ["soc", "nosoc"]

    # Process all data
    print("STEP 1: Reading and processing all data...")
    results = process_all_data(strain_dirs, calc_types)

    print("\n" + "="*60)
    print("STEP 2: Generating all plots...")
    print("="*60 + "\n")

    # Generate all plots
    plot_pdos_comparison(results, strain_dirs)
    plot_dos_comparison_nosoc(results, strain_dirs)
    plot_soc_vs_nosoc_dos(results, strain_dirs)
    plot_vbm_cbm_vs_strain(results, strain_dirs)
    plot_bandgap_vs_strain(results, strain_dirs)

    print("\n" + "="*60)
    print(" All Origin-style plots completed!")
    print("="*60)
    print("\nGenerated files (all at 600 dpi):")
    print("  1. origin_PDOS_comparison.png")
    print("  2. origin_DOS_comparison_noSOC.png")
    print("  3. origin_DOS_SOC_vs_noSOC_strain099.png")
    print("  4. origin_DOS_SOC_vs_noSOC_strain100.png")
    print("  5. origin_DOS_SOC_vs_noSOC_strain101.png")
    print("  6. origin_VBM_CBM_vs_strain.png")
    print("  7. origin_bandgap_vs_strain.png")
    print("\nAll high-resolution images are ready for publication/presentation!")
    print("="*60 + "\n")
