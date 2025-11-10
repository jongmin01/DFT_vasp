#!/usr/bin/env python3
"""
compare_all_dos.py - Compare DOS for all 4 structures
"""

import numpy as np
import matplotlib.pyplot as plt
import os
import sys

def read_doscar(filename='DOSCAR'):
    """Read VASP DOSCAR file"""
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    # Parse header
    header = lines[5].split()
    nedos = int(header[2])
    efermi = float(header[3])
    
    # Read DOS data
    dos_data = []
    for i in range(6, 6 + nedos):
        dos_data.append([float(x) for x in lines[i].split()])
    
    dos_data = np.array(dos_data)
    energy = dos_data[:, 0] - efermi  # Shift to Fermi = 0
    dos_total = dos_data[:, 1]
    
    return energy, dos_total, efermi

def estimate_band_gap(energy, dos, threshold=0.01):
    """Estimate band gap from DOS"""
    below_fermi = energy < 0
    dos_below = dos[below_fermi]
    energy_below = energy[below_fermi]
    
    if len(dos_below) > 0 and np.any(dos_below > threshold):
        vbm = energy_below[dos_below > threshold][-1]
    else:
        vbm = energy_below[-1] if len(energy_below) > 0 else energy[0]
    
    above_fermi = energy > 0
    dos_above = dos[above_fermi]
    energy_above = energy[above_fermi]
    
    if len(dos_above) > 0 and np.any(dos_above > threshold):
        cbm = energy_above[dos_above > threshold][0]
    else:
        cbm = energy_above[0] if len(energy_above) > 0 else energy[-1]
    
    gap = cbm - vbm
    return gap, vbm, cbm

def main():
    """Compare DOS for all structures"""
    
    # Define structures
    structures = [
        'MoSSe_Se_up__WSSe_Se_up',
        'MoSSe_Se_up__WSSe_S_up',
        'MoSSe_S_up__WSSe_Se_up',
        'MoSSe_S_up__WSSe_S_up'
    ]
    
    labels = [
        'Se-up / Se-up',
        'Se-up / S-up',
        'S-up / Se-up',
        'S-up / S-up'
    ]
    
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']
    
    # Read all DOS
    all_data = []
    results = []
    
    print("="*70)
    print("DOS Comparison for MoSSe-WSSe Heterostructures")
    print("="*70)
    print()
    
    for i, struct in enumerate(structures):
        doscar_path = f'{struct}/scf/DOSCAR'
        
        if not os.path.exists(doscar_path):
            print(f"Warning: {doscar_path} not found, skipping...")
            continue
        
        energy, dos, efermi = read_doscar(doscar_path)
        gap, vbm, cbm = estimate_band_gap(energy, dos)
        
        all_data.append((energy, dos, labels[i], colors[i]))
        results.append({
            'structure': struct,
            'label': labels[i],
            'gap': gap,
            'vbm': vbm,
            'cbm': cbm,
            'efermi': efermi
        })
        
        print(f"{labels[i]:20s} | Gap: {gap:5.2f} eV | E_F: {efermi:7.3f} eV")
    
    if len(all_data) == 0:
        print("ERROR: No DOSCAR files found!")
        print("Make sure you're in the parent directory with MoSSe_* subdirectories")
        sys.exit(1)
    
    print("="*70)
    print()
    
    # Plot 1: Individual DOS (2x2 grid)
    fig1, axes = plt.subplots(2, 2, figsize=(14, 10))
    axes = axes.flatten()
    
    for idx, (energy, dos, label, color) in enumerate(all_data):
        ax = axes[idx]
        
        # Plot DOS
        ax.plot(energy, dos, color=color, linewidth=2)
        ax.fill_between(energy, 0, dos, alpha=0.3, color=color)
        
        # Fermi level
        ax.axvline(0, color='k', linestyle='--', linewidth=1.5, alpha=0.7)
        
        # Band gap
        gap = results[idx]['gap']
        vbm = results[idx]['vbm']
        cbm = results[idx]['cbm']
        
        if 0.1 < gap < 10:
            ax.axvline(vbm, color='orange', linestyle=':', linewidth=1.5, alpha=0.7)
            ax.axvline(cbm, color='orange', linestyle=':', linewidth=1.5, alpha=0.7)
            ax.axvspan(vbm, cbm, alpha=0.1, color='gray')
            ax.text(0.5, 0.95, f'Gap ≈ {gap:.2f} eV', 
                   transform=ax.transAxes, ha='center', va='top',
                   fontsize=11, fontweight='bold',
                   bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
        
        # Labels
        ax.set_xlabel('Energy - E$_F$ (eV)', fontsize=12, fontweight='bold')
        ax.set_ylabel('DOS (states/eV)', fontsize=12, fontweight='bold')
        ax.set_title(label, fontsize=13, fontweight='bold')
        ax.set_xlim(-5, 5)
        ax.set_ylim(bottom=0)
        ax.grid(True, alpha=0.3, linestyle='--')
        ax.tick_params(labelsize=10)
    
    plt.tight_layout()
    plt.savefig('DOS_all_structures.png', dpi=300, bbox_inches='tight')
    print("Individual DOS plots saved to: DOS_all_structures.png")
    
    # Plot 2: Overlay comparison
    fig2, ax = plt.subplots(figsize=(10, 7))
    
    for energy, dos, label, color in all_data:
        ax.plot(energy, dos, color=color, linewidth=2, label=label, alpha=0.8)
    
    # Fermi level
    ax.axvline(0, color='k', linestyle='--', linewidth=2, label='E$_F$')
    
    # Labels
    ax.set_xlabel('Energy - E$_F$ (eV)', fontsize=14, fontweight='bold')
    ax.set_ylabel('DOS (states/eV)', fontsize=14, fontweight='bold')
    ax.set_title('DOS Comparison - All Stackings', fontsize=16, fontweight='bold')
    ax.set_xlim(-5, 5)
    ax.set_ylim(bottom=0)
    ax.legend(fontsize=11, framealpha=0.9, loc='best')
    ax.grid(True, alpha=0.3, linestyle='--')
    ax.tick_params(labelsize=12)
    
    plt.tight_layout()
    plt.savefig('DOS_comparison.png', dpi=300, bbox_inches='tight')
    print("Overlay comparison saved to: DOS_comparison.png")
    
    # Save results to file
    with open('DOS_analysis_summary.txt', 'w') as f:
        f.write("="*70 + "\n")
        f.write("DOS Analysis Summary - MoSSe-WSSe Heterostructures\n")
        f.write("="*70 + "\n\n")
        
        f.write(f"{'Structure':<25} {'Band Gap (eV)':<15} {'VBM (eV)':<12} {'CBM (eV)':<12}\n")
        f.write("-"*70 + "\n")
        
        for r in results:
            f.write(f"{r['label']:<25} {r['gap']:>8.3f}       "
                   f"{r['vbm']:>8.3f}    {r['cbm']:>8.3f}\n")
        
        f.write("\n" + "="*70 + "\n")
        f.write("Notes:\n")
        f.write("- Band gaps are estimated from DOS (may differ from actual band structure)\n")
        f.write("- VBM/CBM positions are relative to Fermi level (E_F = 0)\n")
        f.write("- For accurate gap determination, band structure calculation recommended\n")
    
    print("Summary saved to: DOS_analysis_summary.txt")
    print()
    print("Analysis complete! Check the generated files:")
    print("  1. DOS_all_structures.png - Individual plots")
    print("  2. DOS_comparison.png - Overlay comparison")
    print("  3. DOS_analysis_summary.txt - Numerical results")

if __name__ == '__main__':
    main()
