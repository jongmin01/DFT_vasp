#!/usr/bin/env python3
"""
Publication-Quality Multi-Panel Fat Band Structure
Creates a 5-panel figure showing Total + Mo d + W d + S p + Se p contributions

Usage:
    python plot_fatband_multipanel.py

Requirements:
    - PROCAR file in the same directory
    - PyProcar installed (pip install pyprocar)
    - OUTCAR (for Fermi energy extraction)

Author: Generated for MoSSe/WSSe heterostructure analysis
"""

import os
import sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

# Check if PyProcar is available
try:
    import pyprocar
    from pyprocar.core import ElectronicBandStructure
    print("✓ PyProcar loaded successfully")
except ImportError:
    print("✗ Error: PyProcar not found")
    print("  Please install: pip install pyprocar")
    sys.exit(1)

# ============================================================================
# CONFIGURATION
# ============================================================================

# Fermi energy (eV) - will be auto-detected from OUTCAR if available
FERMI_ENERGY = None  # Set to None for auto-detection

# Energy range relative to Fermi level (eV)
ENERGY_MIN = -3.0
ENERGY_MAX = 3.0

# High-symmetry k-points for Γ-M-K-Γ path
KPATH_LABELS = ['$\\Gamma$', 'M', 'K', '$\\Gamma$']
KPATH_TICKS = [0, 43, 83, 120]

# Atom indices (0-indexed) - adjust based on your POSCAR
# Current setup: Mo(0-24), W(25-49), S(50-99), Se(100-149)
ATOMS_MO = list(range(0, 25))
ATOMS_W = list(range(25, 50))
ATOMS_S = list(range(50, 100))
ATOMS_SE = list(range(100, 150))

# d-orbitals: dxy, dyz, dz2, dxz, dx2-y2 (indices 4-8 in VASP)
ORBITALS_D = [4, 5, 6, 7, 8]
# p-orbitals: py, pz, px (indices 1-3 in VASP)
ORBITALS_P = [1, 2, 3]

# Color schemes for better visibility
COLORS = {
    'Mo': 'Reds',
    'W': 'Blues', 
    'S': 'Greens',
    'Se': 'Purples'
}

# Figure settings
DPI = 600
FIGURE_WIDTH = 24  # inches (5 panels × ~4.8 inches each)
FIGURE_HEIGHT = 6
FONT_SIZE = 12
LINE_WIDTH = 1.2

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def extract_fermi_from_outcar(outcar_path='OUTCAR'):
    """
    Extract Fermi energy from OUTCAR file
    
    Returns:
        float: Fermi energy in eV, or None if not found
    """
    if not os.path.exists(outcar_path):
        print(f"  Warning: {outcar_path} not found")
        return None
    
    try:
        with open(outcar_path, 'r') as f:
            for line in f:
                if 'E-fermi' in line:
                    # Format: "E-fermi :   -0.3024     XC(G=0):  -10.8234     alpha+bet : -11.4523"
                    fermi = float(line.split()[2])
                    print(f"  ✓ Fermi energy extracted from OUTCAR: {fermi:.4f} eV")
                    return fermi
    except Exception as e:
        print(f"  Warning: Failed to extract Fermi energy: {e}")
    
    return None

def setup_matplotlib_style():
    """
    Configure matplotlib for publication-quality figures
    """
    plt.rcParams.update({
        'font.size': FONT_SIZE,
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'DejaVu Sans'],
        'axes.linewidth': 1.5,
        'xtick.major.width': 1.5,
        'ytick.major.width': 1.5,
        'xtick.major.size': 6,
        'ytick.major.size': 6,
        'xtick.direction': 'in',
        'ytick.direction': 'in',
        'legend.frameon': False,
        'figure.dpi': 100,
        'savefig.dpi': DPI,
        'savefig.bbox': 'tight',
    })

def check_procar_exists():
    """
    Check if PROCAR file exists
    """
    if os.path.exists('PROCAR'):
        print("✓ PROCAR file found")
        return True
    elif os.path.exists('PROCAR.gz'):
        print("✓ PROCAR.gz file found (will be decompressed by PyProcar)")
        return True
    else:
        print("✗ Error: PROCAR file not found")
        print("  Please make sure PROCAR is in the current directory")
        return False

# ============================================================================
# MAIN PLOTTING FUNCTION
# ============================================================================

def create_multipanel_fatband(fermi=None):
    """
    Create 5-panel fat band structure figure
    
    Args:
        fermi: Fermi energy in eV (if None, will try to extract from OUTCAR)
    """
    # Setup
    setup_matplotlib_style()
    
    # Get Fermi energy
    if fermi is None:
        fermi = extract_fermi_from_outcar()
        if fermi is None:
            print("\n✗ Error: Cannot determine Fermi energy")
            print("  Please either:")
            print("    1. Provide OUTCAR file in current directory")
            print("    2. Set FERMI_ENERGY variable in the script")
            sys.exit(1)
    
    print(f"\nUsing Fermi energy: {fermi:.4f} eV")
    print(f"Energy window: [{ENERGY_MIN}, {ENERGY_MAX}] eV (relative to E_F)")
    
    # Create figure with 5 subplots
    print("\n" + "="*60)
    print("Creating multi-panel fat band structure")
    print("="*60)
    
    fig, axes = plt.subplots(1, 5, figsize=(FIGURE_WIDTH, FIGURE_HEIGHT), 
                             sharey=True, sharex=True)
    
    # Panel configurations
    panel_configs = [
        {
            'ax': axes[0],
            'title': 'Total',
            'mode': 'plain',
            'atoms': None,
            'orbitals': None,
            'cmap': 'viridis',
            'show_colorbar': False
        },
        {
            'ax': axes[1],
            'title': 'Mo d',
            'mode': 'parametric',
            'atoms': ATOMS_MO,
            'orbitals': ORBITALS_D,
            'cmap': COLORS['Mo'],
            'show_colorbar': True
        },
        {
            'ax': axes[2],
            'title': 'W d',
            'mode': 'parametric',
            'atoms': ATOMS_W,
            'orbitals': ORBITALS_D,
            'cmap': COLORS['W'],
            'show_colorbar': True
        },
        {
            'ax': axes[3],
            'title': 'S p',
            'mode': 'parametric',
            'atoms': ATOMS_S,
            'orbitals': ORBITALS_P,
            'cmap': COLORS['S'],
            'show_colorbar': True
        },
        {
            'ax': axes[4],
            'title': 'Se p',
            'mode': 'parametric',
            'atoms': ATOMS_SE,
            'orbitals': ORBITALS_P,
            'cmap': COLORS['Se'],
            'show_colorbar': True
        }
    ]
    
    # Plot each panel
    for i, config in enumerate(panel_configs):
        print(f"\n  Panel {i+1}/5: {config['title']}")
        
        try:
            if config['mode'] == 'plain':
                # Total band structure (no projection)
                pyprocar.bandsplot(
                    code='vasp',
                    dirname='.',
                    mode='plain',
                    fermi=fermi,
                    elimit=[ENERGY_MIN, ENERGY_MAX],
                    kticks=KPATH_TICKS,
                    knames=KPATH_LABELS if i == 0 else None,
                    ax=config['ax'],
                    show=False,
                    linewidth=[LINE_WIDTH, LINE_WIDTH],
                    color=['blue', 'red']
                )
            else:
                # Fat band with projections
                pyprocar.bandsplot(
                    code='vasp',
                    dirname='.',
                    mode='parametric',
                    atoms=config['atoms'],
                    orbitals=config['orbitals'],
                    fermi=fermi,
                    elimit=[ENERGY_MIN, ENERGY_MAX],
                    kticks=KPATH_TICKS,
                    knames=KPATH_LABELS if i == 0 else None,
                    cmap=config['cmap'],
                    vmin=0.0,
                    vmax=1.0,
                    ax=config['ax'],
                    show=False,
                    plot_color_bar=False,  # We'll add custom colorbar
                    linewidth=[LINE_WIDTH, LINE_WIDTH]
                )
                
                # Add custom colorbar for each projection panel
                if config['show_colorbar']:
                    norm = mpl.colors.Normalize(vmin=0, vmax=1)
                    sm = plt.cm.ScalarMappable(cmap=config['cmap'], norm=norm)
                    sm.set_array([])
                    
                    # Position colorbar below each panel
                    cbar_ax = fig.add_axes([
                        axes[i].get_position().x0,
                        0.08,  # Below the main plot
                        axes[i].get_position().width,
                        0.02   # Thin horizontal colorbar
                    ])
                    cbar = plt.colorbar(sm, cax=cbar_ax, orientation='horizontal')
                    cbar.set_label('Weight', fontsize=FONT_SIZE-2)
                    cbar.ax.tick_params(labelsize=FONT_SIZE-3)
            
            # Set title
            config['ax'].set_title(config['title'], fontsize=FONT_SIZE+2, 
                                  fontweight='bold', pad=10)
            
            # Draw Fermi level
            config['ax'].axhline(y=0, color='black', linestyle='--', 
                               linewidth=1.0, alpha=0.7, zorder=1)
            
            # Only show y-label on first panel
            if i == 0:
                config['ax'].set_ylabel('$E - E_F$ (eV)', fontsize=FONT_SIZE+1)
            
            # Grid
            config['ax'].grid(True, alpha=0.3, linestyle=':', linewidth=0.5)
            
            print(f"    ✓ Panel {i+1} complete")
            
        except Exception as e:
            print(f"    ✗ Error in panel {i+1}: {e}")
            continue
    
    # Adjust layout
    plt.tight_layout(rect=[0, 0.12, 1, 0.98])  # Leave space for colorbars
    
    # Save figure
    output_files = [
        'fatband_multipanel.png',
        'fatband_multipanel.pdf',
        'fatband_multipanel.svg'
    ]
    
    print("\n" + "="*60)
    print("Saving figures...")
    print("="*60)
    
    for output_file in output_files:
        try:
            plt.savefig(output_file, dpi=DPI, bbox_inches='tight')
            file_size = os.path.getsize(output_file) / (1024*1024)  # MB
            print(f"  ✓ Saved: {output_file} ({file_size:.2f} MB)")
        except Exception as e:
            print(f"  ✗ Failed to save {output_file}: {e}")
    
    plt.close()
    
    print("\n" + "="*60)
    print("Multi-panel fat band structure generation complete!")
    print("="*60)

# ============================================================================
# ADDITIONAL UTILITY: Individual High-Res Panels
# ============================================================================

def create_individual_panels(fermi=None):
    """
    Create individual high-resolution panels for detailed analysis
    """
    if fermi is None:
        fermi = extract_fermi_from_outcar()
        if fermi is None:
            fermi = FERMI_ENERGY
    
    print("\n" + "="*60)
    print("Creating individual high-resolution panels")
    print("="*60)
    
    individual_configs = [
        ('Total', 'plain', None, None, 'viridis'),
        ('Mo_d', 'parametric', ATOMS_MO, ORBITALS_D, COLORS['Mo']),
        ('W_d', 'parametric', ATOMS_W, ORBITALS_D, COLORS['W']),
        ('S_p', 'parametric', ATOMS_S, ORBITALS_P, COLORS['S']),
        ('Se_p', 'parametric', ATOMS_SE, ORBITALS_P, COLORS['Se']),
    ]
    
    for name, mode, atoms, orbitals, cmap in individual_configs:
        print(f"\n  Creating {name} panel...")
        
        try:
            if mode == 'plain':
                pyprocar.bandsplot(
                    code='vasp',
                    dirname='.',
                    mode='plain',
                    fermi=fermi,
                    elimit=[ENERGY_MIN, ENERGY_MAX],
                    kticks=KPATH_TICKS,
                    knames=KPATH_LABELS,
                    savefig=f'fatband_individual_{name}.png',
                    title=f'Band Structure - {name.replace("_", " ")}',
                    figure_size=(8, 6),
                    dpi=DPI,
                    linewidth=[LINE_WIDTH, LINE_WIDTH]
                )
            else:
                pyprocar.bandsplot(
                    code='vasp',
                    dirname='.',
                    mode='parametric',
                    atoms=atoms,
                    orbitals=orbitals,
                    fermi=fermi,
                    elimit=[ENERGY_MIN, ENERGY_MAX],
                    kticks=KPATH_TICKS,
                    knames=KPATH_LABELS,
                    cmap=cmap,
                    vmin=0.0,
                    vmax=1.0,
                    savefig=f'fatband_individual_{name}.png',
                    title=f'Band Structure - {name.replace("_", " ")}',
                    figure_size=(8, 6),
                    dpi=DPI,
                    linewidth=[LINE_WIDTH, LINE_WIDTH],
                    plot_color_bar=True,
                    colorbar_title='Projection weight'
                )
            
            print(f"    ✓ Saved: fatband_individual_{name}.png")
            
        except Exception as e:
            print(f"    ✗ Error: {e}")

# ============================================================================
# MAIN EXECUTION
# ============================================================================

if __name__ == '__main__':
    print("\n" + "#"*60)
    print("# Multi-Panel Fat Band Structure Generator")
    print("# MoSSe/WSSe Heterostructure Analysis")
    print("#"*60)
    
    # Check if PROCAR exists
    if not check_procar_exists():
        sys.exit(1)
    
    # Create multi-panel figure
    create_multipanel_fatband(fermi=FERMI_ENERGY)
    
    # Optional: Create individual high-res panels
    print("\n" + "="*60)
    user_input = input("Create individual high-resolution panels? (y/n): ").lower()
    if user_input in ['y', 'yes']:
        create_individual_panels(fermi=FERMI_ENERGY)
    
    print("\n" + "#"*60)
    print("# All done! Check your output files.")
    print("#"*60)
    print("\nGenerated files:")
    print("  Main multi-panel:")
    print("    - fatband_multipanel.png (raster, high-res)")
    print("    - fatband_multipanel.pdf (vector, publication)")
    print("    - fatband_multipanel.svg (vector, editable)")
    if user_input in ['y', 'yes']:
        print("\n  Individual panels:")
        print("    - fatband_individual_Total.png")
        print("    - fatband_individual_Mo_d.png")
        print("    - fatband_individual_W_d.png")
        print("    - fatband_individual_S_p.png")
        print("    - fatband_individual_Se_p.png")
    print("\n" + "#"*60)
