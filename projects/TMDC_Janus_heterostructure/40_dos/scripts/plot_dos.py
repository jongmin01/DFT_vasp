#!/usr/bin/env python3
"""
plot_dos.py - Plot DOS from VASP DOSCAR file
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
import os

def read_doscar(filename='DOSCAR'):
    """
    Read VASP DOSCAR file
    
    Returns:
        energy: Energy values (eV)
        dos_total: Total DOS
        dos_up: Spin-up DOS (if spin-polarized)
        dos_down: Spin-down DOS (if spin-polarized)
        efermi: Fermi energy
    """
    with open(filename, 'r') as f:
        # Read header
        lines = f.readlines()
        
    # First line: number of ions, etc.
    # Line 6: NEDOS, Fermi energy
    header = lines[5].split()
    nedos = int(header[2])
    efermi = float(header[3])
    
    # Read DOS data (starts at line 6)
    dos_data = []
    for i in range(6, 6 + nedos):
        dos_data.append([float(x) for x in lines[i].split()])
    
    dos_data = np.array(dos_data)
    
    energy = dos_data[:, 0] - efermi  # Shift to Fermi level = 0
    
    # Check if spin-polarized
    if dos_data.shape[1] == 3:  # Non-spin-polarized
        dos_total = dos_data[:, 1]
        dos_up = None
        dos_down = None
    elif dos_data.shape[1] == 5:  # Spin-polarized
        dos_up = dos_data[:, 1]
        dos_down = -dos_data[:, 2]  # Negative for spin-down
        dos_total = dos_data[:, 1] + dos_data[:, 2]
    else:
        raise ValueError("Unexpected DOSCAR format")
    
    return energy, dos_total, dos_up, dos_down, efermi

def estimate_band_gap(energy, dos, threshold=0.01):
    """
    Estimate band gap from DOS
    
    Args:
        energy: Energy array
        dos: DOS array
        threshold: DOS threshold for gap detection
    
    Returns:
        gap: Band gap (eV)
        vbm: Valence band maximum
        cbm: Conduction band minimum
    """
    # Find indices near Fermi level (E=0)
    fermi_idx = np.argmin(np.abs(energy))
    
    # Find VBM (highest energy with DOS > threshold below Fermi)
    below_fermi = energy < 0
    dos_below = dos[below_fermi]
    energy_below = energy[below_fermi]
    
    if len(dos_below) > 0:
        occupied = dos_below > threshold
        if np.any(occupied):
            vbm = energy_below[occupied][-1]
        else:
            vbm = energy_below[-1]
    else:
        vbm = energy[0]
    
    # Find CBM (lowest energy with DOS > threshold above Fermi)
    above_fermi = energy > 0
    dos_above = dos[above_fermi]
    energy_above = energy[above_fermi]
    
    if len(dos_above) > 0:
        unoccupied = dos_above > threshold
        if np.any(unoccupied):
            cbm = energy_above[unoccupied][0]
        else:
            cbm = energy_above[0]
    else:
        cbm = energy[-1]
    
    gap = cbm - vbm
    
    return gap, vbm, cbm

def plot_dos(energy, dos_total, dos_up=None, dos_down=None, 
             title='Density of States', output='DOS.png'):
    """
    Plot DOS
    
    Args:
        energy: Energy array (shifted to Fermi level = 0)
        dos_total: Total DOS
        dos_up: Spin-up DOS (optional)
        dos_down: Spin-down DOS (optional)
        title: Plot title
        output: Output filename
    """
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Plot DOS
    if dos_up is not None and dos_down is not None:
        # Spin-polarized
        ax.plot(energy, dos_up, 'r-', linewidth=2, label='Spin up')
        ax.plot(energy, dos_down, 'b-', linewidth=2, label='Spin down')
        ax.fill_between(energy, 0, dos_up, alpha=0.3, color='r')
        ax.fill_between(energy, 0, dos_down, alpha=0.3, color='b')
    else:
        # Non-spin or total
        ax.plot(energy, dos_total, 'b-', linewidth=2, label='Total DOS')
        ax.fill_between(energy, 0, dos_total, alpha=0.3, color='b')
    
    # Fermi level
    ax.axvline(0, color='k', linestyle='--', linewidth=1.5, label='Fermi level')
    
    # Estimate and show band gap
    gap, vbm, cbm = estimate_band_gap(energy, dos_total)
    if 0.1 < gap < 10:  # Reasonable gap range
        ax.axvline(vbm, color='orange', linestyle=':', linewidth=1.5, alpha=0.7)
        ax.axvline(cbm, color='orange', linestyle=':', linewidth=1.5, alpha=0.7)
        ax.text(0.02, 0.98, f'Band gap ≈ {gap:.2f} eV', 
                transform=ax.transAxes, fontsize=12,
                verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    # Labels and formatting
    ax.set_xlabel('Energy - E$_F$ (eV)', fontsize=14, fontweight='bold')
    ax.set_ylabel('DOS (states/eV)', fontsize=14, fontweight='bold')
    ax.set_title(title, fontsize=16, fontweight='bold')
    ax.set_xlim(-5, 5)
    ax.set_ylim(bottom=0)
    ax.legend(fontsize=12, framealpha=0.9)
    ax.grid(True, alpha=0.3, linestyle='--')
    ax.tick_params(labelsize=12)
    
    plt.tight_layout()
    plt.savefig(output, dpi=300, bbox_inches='tight')
    print(f"DOS plot saved to: {output}")
    
    return gap, vbm, cbm

def main():
    """Main function"""
    # Check if DOSCAR exists
    if len(sys.argv) > 1:
        doscar_file = sys.argv[1]
    else:
        doscar_file = 'DOSCAR'
    
    if not os.path.exists(doscar_file):
        print(f"ERROR: {doscar_file} not found!")
        print("Usage: python plot_dos.py [DOSCAR]")
        sys.exit(1)
    
    # Get directory name for title
    dirname = os.path.basename(os.path.dirname(os.path.abspath(doscar_file)))
    if dirname == '':
        dirname = 'DOS'
    
    print(f"Reading DOS from: {doscar_file}")
    
    # Read DOSCAR
    try:
        energy, dos_total, dos_up, dos_down, efermi = read_doscar(doscar_file)
        print(f"  NEDOS: {len(energy)}")
        print(f"  Fermi energy (original): {efermi:.4f} eV")
        print(f"  Energy range: {energy[0]:.2f} to {energy[-1]:.2f} eV (shifted)")
    except Exception as e:
        print(f"ERROR reading DOSCAR: {e}")
        sys.exit(1)
    
    # Plot
    output_file = doscar_file.replace('DOSCAR', 'DOS.png')
    if output_file == doscar_file:
        output_file = 'DOS.png'
    
    gap, vbm, cbm = plot_dos(energy, dos_total, dos_up, dos_down, 
                              title=f'DOS - {dirname}', output=output_file)
    
    # Print results
    print("\n" + "="*50)
    print("DOS Analysis Results")
    print("="*50)
    print(f"Band gap (estimated): {gap:.3f} eV")
    print(f"VBM (relative to E_F): {vbm:.3f} eV")
    print(f"CBM (relative to E_F): {cbm:.3f} eV")
    
    if gap < 0.1:
        print("\nNote: Very small or no band gap detected")
        print("      This might be a metal or the gap is too small")
    elif gap > 5:
        print("\nNote: Very large band gap detected")
        print("      Check if this is reasonable for your system")
    
    print("="*50)

if __name__ == '__main__':
    main()
