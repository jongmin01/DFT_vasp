#!/usr/bin/env python3
"""
plot_band_structure.py - Plot band structure from VASP EIGENVAL
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
import os

def read_eigenval(filename='EIGENVAL'):
    """
    Read VASP EIGENVAL file
    
    Returns:
        kpoints: k-point coordinates
        energies: eigenvalues [nkpts, nbands, nspin]
        efermi: Fermi energy
    """
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    # Skip header (5 lines)
    # Line 6: nelect, nkpts, nbands
    header = lines[5].split()
    nkpts = int(header[1])
    nbands = int(header[2])
    
    # Determine if spin-polarized
    # Check line 7 format
    test_line = lines[7].split()
    if len(test_line) == 4:  # k-point line: kx ky kz weight
        data_start = 7
        spin_polarized = False
    else:
        data_start = 8
        spin_polarized = False
    
    # Read data
    kpoints = []
    energies = []
    
    i = data_start
    for _ in range(nkpts):
        # k-point line
        kpt_line = lines[i].split()
        kpoints.append([float(kpt_line[0]), float(kpt_line[1]), float(kpt_line[2])])
        i += 1
        
        # Band energies
        bands = []
        for _ in range(nbands):
            band_line = lines[i].split()
            # Format: band_index energy [occ] [energy2 occ2] for spin
            bands.append(float(band_line[1]))
            i += 1
        
        energies.append(bands)
        i += 1  # Skip blank line
    
    kpoints = np.array(kpoints)
    energies = np.array(energies)  # [nkpts, nbands]
    
    return kpoints, energies

def read_fermi_from_doscar(doscar='DOSCAR'):
    """Read Fermi energy from DOSCAR"""
    try:
        with open(doscar, 'r') as f:
            lines = f.readlines()
        header = lines[5].split()
        efermi = float(header[3])
        return efermi
    except:
        return 0.0

def get_kpath_distance(kpoints):
    """Calculate cumulative distance along k-path"""
    distances = [0.0]
    for i in range(1, len(kpoints)):
        dk = np.linalg.norm(kpoints[i] - kpoints[i-1])
        distances.append(distances[-1] + dk)
    return np.array(distances)

def find_band_gap(energies, efermi=0.0):
    """
    Find band gap and determine if direct or indirect
    
    Returns:
        gap: band gap (eV)
        vbm: valence band maximum
        cbm: conduction band minimum
        vbm_k: k-point index of VBM
        cbm_k: k-point index of CBM
        is_direct: True if direct gap
    """
    # Shift energies to Fermi level
    energies_shifted = energies - efermi
    
    # Find VBM (highest occupied state)
    occupied = energies_shifted < 0
    if not np.any(occupied):
        return None, None, None, None, None, None
    
    vbm_indices = np.where(occupied)
    vbm_band = vbm_indices[1].max()
    vbm_k = vbm_indices[0][vbm_indices[1] == vbm_band][0]
    vbm = energies_shifted[vbm_k, vbm_band]
    
    # Find CBM (lowest unoccupied state)
    unoccupied = energies_shifted > 0
    if not np.any(unoccupied):
        return None, None, None, None, None, None
    
    cbm_indices = np.where(unoccupied)
    cbm_band = cbm_indices[1].min()
    cbm_k = cbm_indices[0][cbm_indices[1] == cbm_band][0]
    cbm = energies_shifted[cbm_k, cbm_band]
    
    gap = cbm - vbm
    is_direct = (vbm_k == cbm_k)
    
    return gap, vbm, cbm, vbm_k, cbm_k, is_direct

def plot_band_structure(kpoints, energies, efermi=0.0, 
                       high_symmetry_points=None,
                       high_symmetry_labels=None,
                       title='Band Structure',
                       output='band_structure.png'):
    """
    Plot band structure
    
    Args:
        kpoints: k-point coordinates
        energies: eigenvalues [nkpts, nbands]
        efermi: Fermi energy
        high_symmetry_points: list of k-point indices for high-symmetry points
        high_symmetry_labels: list of labels for high-symmetry points
        title: plot title
        output: output filename
    """
    # Get k-path distance
    kdist = get_kpath_distance(kpoints)
    
    # Shift energies to Fermi level
    energies_shifted = energies - efermi
    
    # Find band gap
    gap, vbm, cbm, vbm_k, cbm_k, is_direct = find_band_gap(energies, efermi)
    
    # Plot
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Plot all bands
    for iband in range(energies.shape[1]):
        ax.plot(kdist, energies_shifted[:, iband], 'b-', linewidth=1.0, alpha=0.7)
    
    # Fermi level
    ax.axhline(0, color='k', linestyle='--', linewidth=1.5, alpha=0.7, label='E$_F$')
    
    # Mark VBM and CBM if gap exists
    if gap is not None and 0.1 < gap < 5:
        ax.plot(kdist[vbm_k], vbm, 'ro', markersize=8, label=f'VBM')
        ax.plot(kdist[cbm_k], cbm, 'go', markersize=8, label=f'CBM')
        
        gap_type = "Direct" if is_direct else "Indirect"
        ax.text(0.02, 0.98, f'Band gap: {gap:.3f} eV ({gap_type})', 
                transform=ax.transAxes, fontsize=12,
                verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    # High-symmetry points
    if high_symmetry_points is not None and high_symmetry_labels is not None:
        for kpt, label in zip(high_symmetry_points, high_symmetry_labels):
            ax.axvline(kdist[kpt], color='gray', linestyle='-', linewidth=0.5, alpha=0.5)
        
        ax.set_xticks([kdist[kpt] for kpt in high_symmetry_points])
        ax.set_xticklabels(high_symmetry_labels, fontsize=12)
    
    # Labels and formatting
    ax.set_ylabel('Energy - E$_F$ (eV)', fontsize=14, fontweight='bold')
    ax.set_title(title, fontsize=16, fontweight='bold')
    ax.set_xlim(kdist[0], kdist[-1])
    ax.set_ylim(-3, 3)
    ax.legend(fontsize=10, loc='upper right')
    ax.grid(True, alpha=0.3, linestyle='--', axis='y')
    ax.tick_params(labelsize=12)
    
    plt.tight_layout()
    plt.savefig(output, dpi=300, bbox_inches='tight')
    print(f"Band structure plot saved to: {output}")
    
    return gap, vbm, cbm, vbm_k, cbm_k, is_direct

def main():
    """Main function"""
    # Check if EIGENVAL exists
    if len(sys.argv) > 1:
        eigenval_file = sys.argv[1]
        dirname = os.path.dirname(eigenval_file)
    else:
        eigenval_file = 'EIGENVAL'
        dirname = '.'
    
    if not os.path.exists(eigenval_file):
        print(f"ERROR: {eigenval_file} not found!")
        print("Usage: python plot_band_structure.py [EIGENVAL]")
        sys.exit(1)
    
    print(f"Reading band structure from: {eigenval_file}")
    
    # Read EIGENVAL
    try:
        kpoints, energies = read_eigenval(eigenval_file)
        print(f"  Number of k-points: {len(kpoints)}")
        print(f"  Number of bands: {energies.shape[1]}")
    except Exception as e:
        print(f"ERROR reading EIGENVAL: {e}")
        sys.exit(1)
    
    # Read Fermi energy
    doscar_file = os.path.join(dirname, 'DOSCAR')
    if os.path.exists(doscar_file):
        efermi = read_fermi_from_doscar(doscar_file)
        print(f"  Fermi energy: {efermi:.4f} eV")
    else:
        efermi = 0.0
        print(f"  DOSCAR not found, assuming E_F = 0")
    
    # Define high-symmetry points (Γ-M-K-Γ)
    # Assuming 40 points per segment
    npts_per_segment = 40
    high_sym_points = [0, npts_per_segment, 2*npts_per_segment, 3*npts_per_segment]
    high_sym_labels = ['Γ', 'M', 'K', 'Γ']
    
    # Auto-detect number of segments
    nkpts = len(kpoints)
    if nkpts != 3 * npts_per_segment + 1:
        # Recalculate
        npts_per_segment = nkpts // 3
        high_sym_points = [0, npts_per_segment, 2*npts_per_segment, nkpts-1]
    
    # Get structure name
    struct_name = os.path.basename(os.path.dirname(os.path.dirname(eigenval_file)))
    if struct_name == '':
        struct_name = 'Band Structure'
    
    # Plot
    output_file = eigenval_file.replace('EIGENVAL', 'band_structure.png')
    if output_file == eigenval_file:
        output_file = 'band_structure.png'
    
    gap, vbm, cbm, vbm_k, cbm_k, is_direct = plot_band_structure(
        kpoints, energies, efermi,
        high_symmetry_points=high_sym_points,
        high_symmetry_labels=high_sym_labels,
        title=f'Band Structure - {struct_name} (SOC OFF)',
        output=output_file
    )
    
    # Print results
    print("\n" + "="*60)
    print("Band Structure Analysis")
    print("="*60)
    
    if gap is not None:
        print(f"Band gap: {gap:.3f} eV")
        print(f"VBM: {vbm:.3f} eV (k-point {vbm_k})")
        print(f"CBM: {cbm:.3f} eV (k-point {cbm_k})")
        print(f"Gap type: {'Direct' if is_direct else 'Indirect'}")
        
        if is_direct:
            print(f"  → Direct gap at k-point {vbm_k}")
        else:
            print(f"  → Indirect gap:")
            print(f"    VBM at k-point {vbm_k}")
            print(f"    CBM at k-point {cbm_k}")
    else:
        print("Could not determine band gap")
    
    print("="*60)

if __name__ == '__main__':
    main()
