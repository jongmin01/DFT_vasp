#!/usr/bin/env python3
"""
Improved Publication-Quality Multi-Panel Fat Band Structure
- Fixed energy range for all panels
- Clear panel titles and labels
- Enhanced colorbar visibility
- Better visualization of weak contributions

Usage:
    python plot_fatband_multipanel_improved.py

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
from matplotlib.gridspec import GridSpec

# Check if PyProcar is available
try:
    import pyprocar
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

# Color schemes - optimized for better visibility
COLORS = {
    'Mo': 'Reds',      # Red for Mo
    'W': 'Blues',      # Blue for W
    'S': 'Greens',     # Green for S
    'Se': 'Purples'    # Purple for Se
}

# Figure settings
DPI = 600
FIGURE_WIDTH = 25  # inches (5 panels × ~5 inches each)
FIGURE_HEIGHT = 7
FONT_SIZE = 14
TITLE_FONT_SIZE = 16
LINE_WIDTH = 1.5

# Colorbar settings - adjusted for weak contributions
VMAX_METAL = 1.0      # For Mo d, W d (strong contributions)
VMAX_CHALCOGEN = 0.3  # For S p, Se p (weak contributions, amplified)

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
        'font.sans-serif': ['Arial', 'DejaVu Sans', 'Helvetica'],
        'axes.linewidth': 1.8,
        'axes.labelsize': FONT_SIZE,
        'xtick.major.width': 1.5,
        'ytick.major.width': 1.5,
        'xtick.major.size': 6,
        'ytick.major.size': 6,
        'xtick.minor.size': 4,
        'ytick.minor.size': 4,
        'xtick.direction': 'in',
        'ytick.direction': 'in',
        'xtick.top': True,
        'ytick.right': True,
        'legend.frameon': False,
        'figure.dpi': 100,
        'savefig.dpi': DPI,
        'savefig.bbox': 'tight',
    })

def check_files():
    """
    Check if required files exist
    """
    issues = []
    
    if os.path.exists('PROCAR'):
        print("✓ PROCAR file found")
    elif os.path.exists('PROCAR.gz'):
        print("✓ PROCAR.gz file found")
    else:
        issues.append("PROCAR file not found")
    
    if not os.path.exists('OUTCAR'):
        print("⚠ OUTCAR not found - will need manual Fermi energy")
    
    return len(issues) == 0

# ============================================================================
# MAIN PLOTTING FUNCTION - IMPROVED VERSION
# ============================================================================

def create_improved_multipanel_fatband(fermi=None):
    """
    Create improved 5-panel fat band structure figure with better visibility
    
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
    
    # Create figure with GridSpec for better control
    print("\n" + "="*70)
    print("Creating improved multi-panel fat band structure")
    print("="*70)
    
    fig = plt.figure(figsize=(FIGURE_WIDTH, FIGURE_HEIGHT))
    gs = GridSpec(1, 5, figure=fig, wspace=0.08, hspace=0.05,
                  left=0.05, right=0.98, top=0.92, bottom=0.18)
    
    axes = [fig.add_subplot(gs[0, i]) for i in range(5)]
    
    # Panel configurations with improved settings
    panel_configs = [
        {
            'ax': axes[0],
            'title': 'Total',
            'subtitle': 'Band Structure',
            'mode': 'plain',
            'atoms': None,
            'orbitals': None,
            'cmap': 'viridis',
            'vmax': 1.0,
            'show_colorbar': False
        },
        {
            'ax': axes[1],
            'title': 'Mo d orbitals',
            'subtitle': 'MoSSe Layer',
            'mode': 'parametric',
            'atoms': ATOMS_MO,
            'orbitals': ORBITALS_D,
            'cmap': COLORS['Mo'],
            'vmax': VMAX_METAL,
            'show_colorbar': True
        },
        {
            'ax': axes[2],
            'title': 'W d orbitals',
            'subtitle': 'WSSe Layer',
            'mode': 'parametric',
            'atoms': ATOMS_W,
            'orbitals': ORBITALS_D,
            'cmap': COLORS['W'],
            'vmax': VMAX_METAL,
            'show_colorbar': True
        },
        {
            'ax': axes[3],
            'title': 'S p orbitals',
            'subtitle': 'Chalcogen',
            'mode': 'parametric',
            'atoms': ATOMS_S,
            'orbitals': ORBITALS_P,
            'cmap': COLORS['S'],
            'vmax': VMAX_CHALCOGEN,  # Lower vmax to amplify weak contributions
            'show_colorbar': True
        },
        {
            'ax': axes[4],
            'title': 'Se p orbitals',
            'subtitle': 'Chalcogen',
            'mode': 'parametric',
            'atoms': ATOMS_SE,
            'orbitals': ORBITALS_P,
            'cmap': COLORS['Se'],
            'vmax': VMAX_CHALCOGEN,  # Lower vmax to amplify weak contributions
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
                    knames=KPATH_LABELS,
                    ax=config['ax'],
                    show=False,
                    linewidth=[LINE_WIDTH, LINE_WIDTH],
                    color=['#D62728', '#1F77B4']  # Red for spin up, blue for spin down
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
                    knames=KPATH_LABELS,
                    cmap=config['cmap'],
                    vmin=0.0,
                    vmax=config['vmax'],
                    ax=config['ax'],
                    show=False,
                    plot_color_bar=False,  # We'll add custom colorbar
                    linewidth=[LINE_WIDTH, LINE_WIDTH]
                )
            
            # Set energy limits explicitly
            config['ax'].set_ylim([ENERGY_MIN, ENERGY_MAX])
            
            # Set title with two lines for clarity
            title_text = f"{config['title']}\n{config['subtitle']}"
            config['ax'].set_title(title_text, fontsize=TITLE_FONT_SIZE, 
                                  fontweight='bold', pad=12, linespacing=1.3)
            
            # Draw Fermi level - more prominent
            config['ax'].axhline(y=0, color='black', linestyle='--', 
                               linewidth=1.5, alpha=0.8, zorder=10)
            
            # Y-axis settings
            if i == 0:
                config['ax'].set_ylabel('$E - E_\\mathrm{F}$ (eV)', 
                                       fontsize=FONT_SIZE+2, fontweight='bold')
                # Set y-ticks
                yticks = np.arange(ENERGY_MIN, ENERGY_MAX + 0.5, 1.0)
                config['ax'].set_yticks(yticks)
            else:
                config['ax'].set_yticklabels([])
            
            # X-axis label for middle panel
            if i == 2:
                config['ax'].set_xlabel('Wave Vector', fontsize=FONT_SIZE+2, 
                                       fontweight='bold')
            
            # Grid with subtle styling
            config['ax'].grid(True, alpha=0.25, linestyle=':', linewidth=0.8, 
                            color='gray', zorder=0)
            
            # Add colorbar for projection panels
            if config['show_colorbar']:
                # Create colorbar below each panel
                cbar_ax = fig.add_axes([
                    axes[i].get_position().x0 + 0.005,  # Slight offset from left
                    0.08,  # Below the main plot
                    axes[i].get_position().width - 0.01,  # Slightly narrower
                    0.025   # Thin horizontal colorbar
                ])
                
                norm = mpl.colors.Normalize(vmin=0, vmax=config['vmax'])
                sm = plt.cm.ScalarMappable(cmap=config['cmap'], norm=norm)
                sm.set_array([])
                
                cbar = plt.colorbar(sm, cax=cbar_ax, orientation='horizontal')
                
                # Colorbar label
                if config['vmax'] == VMAX_METAL:
                    cbar.set_label('Projection Weight', fontsize=FONT_SIZE-1, 
                                  fontweight='bold')
                else:
                    cbar.set_label('Projection Weight (×3 amplified)', 
                                  fontsize=FONT_SIZE-1, fontweight='bold')
                
                # Colorbar ticks
                cbar.ax.tick_params(labelsize=FONT_SIZE-2)
                if config['vmax'] == VMAX_METAL:
                    cbar.set_ticks([0, 0.25, 0.5, 0.75, 1.0])
                else:
                    cbar.set_ticks([0, 0.1, 0.2, 0.3])
                
                # Add frame to colorbar
                cbar.outline.set_linewidth(1.5)
            
            print(f"    ✓ Panel {i+1} complete")
            
        except Exception as e:
            print(f"    ✗ Error in panel {i+1}: {e}")
            import traceback
            traceback.print_exc()
            continue
    
    # Add overall figure title
    fig.suptitle('MoSSe/WSSe Heterostructure: Orbital-Resolved Band Structure',
                 fontsize=TITLE_FONT_SIZE+2, fontweight='bold', y=0.98)
    
    # Save figures
    output_files = [
        'fatband_multipanel_improved.png',
        'fatband_multipanel_improved.pdf',
        'fatband_multipanel_improved.svg'
    ]
    
    print("\n" + "="*70)
    print("Saving figures...")
    print("="*70)
    
    for output_file in output_files:
        try:
            plt.savefig(output_file, dpi=DPI, bbox_inches='tight', 
                       facecolor='white', edgecolor='none')
            file_size = os.path.getsize(output_file) / (1024*1024)  # MB
            print(f"  ✓ Saved: {output_file} ({file_size:.2f} MB)")
        except Exception as e:
            print(f"  ✗ Failed to save {output_file}: {e}")
    
    plt.close()
    
    print("\n" + "="*70)
    print("Improved multi-panel fat band structure generation complete!")
    print("="*70)

# ============================================================================
# ALTERNATIVE: 3-PANEL VERSION (Only important orbitals)
# ============================================================================

def create_compact_3panel_fatband(fermi=None):
    """
    Create compact 3-panel version focusing on important orbitals
    [Total] [Mo d] [W d]
    
    Omits S/Se p orbitals which have negligible contribution
    """
    if fermi is None:
        fermi = extract_fermi_from_outcar()
        if fermi is None:
            fermi = FERMI_ENERGY
    
    print("\n" + "="*70)
    print("Creating compact 3-panel fat band structure (Total + Mo d + W d)")
    print("="*70)
    
    setup_matplotlib_style()
    
    fig = plt.figure(figsize=(18, 7))
    gs = GridSpec(1, 3, figure=fig, wspace=0.08,
                  left=0.06, right=0.98, top=0.92, bottom=0.18)
    
    axes = [fig.add_subplot(gs[0, i]) for i in range(3)]
    
    panel_configs = [
        {
            'ax': axes[0],
            'title': 'Total Band Structure',
            'mode': 'plain',
            'atoms': None,
            'orbitals': None,
            'cmap': 'viridis',
            'show_colorbar': False
        },
        {
            'ax': axes[1],
            'title': 'Mo d orbitals\n(MoSSe Layer)',
            'mode': 'parametric',
            'atoms': ATOMS_MO,
            'orbitals': ORBITALS_D,
            'cmap': COLORS['Mo'],
            'show_colorbar': True
        },
        {
            'ax': axes[2],
            'title': 'W d orbitals\n(WSSe Layer)',
            'mode': 'parametric',
            'atoms': ATOMS_W,
            'orbitals': ORBITALS_D,
            'cmap': COLORS['W'],
            'show_colorbar': True
        }
    ]
    
    for i, config in enumerate(panel_configs):
        print(f"\n  Panel {i+1}/3: {config['title'].split()[0]}")
        
        try:
            if config['mode'] == 'plain':
                pyprocar.bandsplot(
                    code='vasp',
                    dirname='.',
                    mode='plain',
                    fermi=fermi,
                    elimit=[ENERGY_MIN, ENERGY_MAX],
                    kticks=KPATH_TICKS,
                    knames=KPATH_LABELS,
                    ax=config['ax'],
                    show=False,
                    linewidth=[LINE_WIDTH, LINE_WIDTH],
                    color=['#D62728', '#1F77B4']
                )
            else:
                pyprocar.bandsplot(
                    code='vasp',
                    dirname='.',
                    mode='parametric',
                    atoms=config['atoms'],
                    orbitals=config['orbitals'],
                    fermi=fermi,
                    elimit=[ENERGY_MIN, ENERGY_MAX],
                    kticks=KPATH_TICKS,
                    knames=KPATH_LABELS,
                    cmap=config['cmap'],
                    vmin=0.0,
                    vmax=1.0,
                    ax=config['ax'],
                    show=False,
                    plot_color_bar=False,
                    linewidth=[LINE_WIDTH, LINE_WIDTH]
                )
            
            config['ax'].set_ylim([ENERGY_MIN, ENERGY_MAX])
            config['ax'].set_title(config['title'], fontsize=TITLE_FONT_SIZE+2, 
                                  fontweight='bold', pad=12, linespacing=1.3)
            config['ax'].axhline(y=0, color='black', linestyle='--', 
                               linewidth=1.5, alpha=0.8, zorder=10)
            
            if i == 0:
                config['ax'].set_ylabel('$E - E_\\mathrm{F}$ (eV)', 
                                       fontsize=FONT_SIZE+3, fontweight='bold')
                yticks = np.arange(ENERGY_MIN, ENERGY_MAX + 0.5, 1.0)
                config['ax'].set_yticks(yticks)
            else:
                config['ax'].set_yticklabels([])
            
            if i == 1:
                config['ax'].set_xlabel('Wave Vector', fontsize=FONT_SIZE+3, 
                                       fontweight='bold')
            
            config['ax'].grid(True, alpha=0.25, linestyle=':', linewidth=0.8, 
                            color='gray', zorder=0)
            
            if config['show_colorbar']:
                cbar_ax = fig.add_axes([
                    axes[i].get_position().x0 + 0.01,
                    0.08,
                    axes[i].get_position().width - 0.02,
                    0.03
                ])
                norm = mpl.colors.Normalize(vmin=0, vmax=1.0)
                sm = plt.cm.ScalarMappable(cmap=config['cmap'], norm=norm)
                sm.set_array([])
                cbar = plt.colorbar(sm, cax=cbar_ax, orientation='horizontal')
                cbar.set_label('Projection Weight', fontsize=FONT_SIZE, 
                              fontweight='bold')
                cbar.ax.tick_params(labelsize=FONT_SIZE-1)
                cbar.set_ticks([0, 0.25, 0.5, 0.75, 1.0])
                cbar.outline.set_linewidth(1.5)
            
            print(f"    ✓ Panel {i+1} complete")
            
        except Exception as e:
            print(f"    ✗ Error: {e}")
    
    fig.suptitle('MoSSe/WSSe Heterostructure: d-Orbital Dominated Band Structure',
                 fontsize=TITLE_FONT_SIZE+3, fontweight='bold', y=0.98)
    
    output_files = [
        'fatband_3panel_compact.png',
        'fatband_3panel_compact.pdf'
    ]
    
    print("\n" + "="*70)
    print("Saving compact 3-panel figures...")
    print("="*70)
    
    for output_file in output_files:
        try:
            plt.savefig(output_file, dpi=DPI, bbox_inches='tight',
                       facecolor='white', edgecolor='none')
            file_size = os.path.getsize(output_file) / (1024*1024)
            print(f"  ✓ Saved: {output_file} ({file_size:.2f} MB)")
        except Exception as e:
            print(f"  ✗ Failed: {e}")
    
    plt.close()
    print("\n✓ Compact 3-panel version complete!")

# ============================================================================
# MAIN EXECUTION
# ============================================================================

if __name__ == '__main__':
    print("\n" + "#"*70)
    print("# Improved Multi-Panel Fat Band Structure Generator")
    print("# MoSSe/WSSe Heterostructure Analysis")
    print("#"*70)
    
    # Check if required files exist
    if not check_files():
        sys.exit(1)
    
    # Create improved 5-panel figure
    create_improved_multipanel_fatband(fermi=FERMI_ENERGY)
    
    # Optional: Create compact 3-panel version
    print("\n" + "="*70)
    user_input = input("Create compact 3-panel version (Total + Mo d + W d only)? (y/n): ").lower()
    if user_input in ['y', 'yes']:
        create_compact_3panel_fatband(fermi=FERMI_ENERGY)
    
    print("\n" + "#"*70)
    print("# All done! Check your output files.")
    print("#"*70)
    print("\nGenerated files:")
    print("  Main 5-panel (improved):")
    print("    - fatband_multipanel_improved.png")
    print("    - fatband_multipanel_improved.pdf")
    print("    - fatband_multipanel_improved.svg")
    if user_input in ['y', 'yes']:
        print("\n  Compact 3-panel (d-orbitals focused):")
        print("    - fatband_3panel_compact.png")
        print("    - fatband_3panel_compact.pdf")
    print("\n" + "#"*70)
    print("\nKey improvements:")
    print("  ✓ Fixed energy range [-3, 3] eV for all panels")
    print("  ✓ Clear panel titles with subtitles")
    print("  ✓ Enhanced colorbar visibility and labeling")
    print("  ✓ S/Se p orbitals amplified (×3) for better visibility")
    print("  ✓ More prominent Fermi level line")
    print("  ✓ Better grid and axis styling")
    print("  ✓ Optional compact 3-panel version (main orbitals only)")
    print("#"*70)
