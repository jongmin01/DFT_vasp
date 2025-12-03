#!/usr/bin/env python3
"""
Robust Multi-Panel Fat Band Structure - Alternative Version
Uses separate plotting calls to avoid PyProcar colorbar bugs

Usage:
    python plot_fatband_multipanel_robust.py

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

FERMI_ENERGY = None  # Auto-detect from OUTCAR
ENERGY_MIN = -3.0
ENERGY_MAX = 3.0

KPATH_LABELS = ['$\\Gamma$', 'M', 'K', '$\\Gamma$']
KPATH_TICKS = [0, 43, 83, 120]

ATOMS_MO = list(range(0, 25))
ATOMS_W = list(range(25, 50))
ATOMS_S = list(range(50, 100))
ATOMS_SE = list(range(100, 150))

ORBITALS_D = [4, 5, 6, 7, 8]
ORBITALS_P = [1, 2, 3]

COLORS = {
    'Mo': 'Reds',
    'W': 'Blues', 
    'S': 'Greens',
    'Se': 'Purples'
}

DPI = 600
FIGURE_WIDTH = 25
FIGURE_HEIGHT = 7
FONT_SIZE = 14
TITLE_FONT_SIZE = 16
LINE_WIDTH = 1.5

VMAX_METAL = 1.0
VMAX_CHALCOGEN = 0.3

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def extract_fermi_from_outcar(outcar_path='OUTCAR'):
    if not os.path.exists(outcar_path):
        return None
    try:
        with open(outcar_path, 'r') as f:
            for line in f:
                if 'E-fermi' in line:
                    fermi = float(line.split()[2])
                    print(f"  ✓ Fermi energy: {fermi:.4f} eV")
                    return fermi
    except:
        pass
    return None

def setup_matplotlib():
    plt.rcParams.update({
        'font.size': FONT_SIZE,
        'font.family': 'sans-serif',
        'axes.linewidth': 1.8,
        'xtick.major.width': 1.5,
        'ytick.major.width': 1.5,
        'xtick.major.size': 6,
        'ytick.major.size': 6,
        'xtick.direction': 'in',
        'ytick.direction': 'in',
        'xtick.top': True,
        'ytick.right': True,
        'savefig.dpi': DPI,
        'savefig.bbox': 'tight',
    })

def check_files():
    if os.path.exists('PROCAR') or os.path.exists('PROCAR.gz'):
        print("✓ PROCAR file found")
        return True
    print("✗ PROCAR file not found")
    return False

# ============================================================================
# ROBUST PLOTTING - Each panel plotted separately
# ============================================================================

def plot_single_panel(ax, title, subtitle, mode, atoms, orbitals, cmap, vmax, fermi, show_ylabel=False):
    """
    Plot a single panel with proper error handling
    """
    try:
        if mode == 'plain':
            # Total band - simple mode
            pyprocar.bandsplot(
                code='vasp',
                dirname='.',
                mode='plain',
                fermi=fermi,
                elimit=[ENERGY_MIN, ENERGY_MAX],
                kticks=KPATH_TICKS,
                knames=KPATH_LABELS,
                ax=ax,
                show=False,
                savefig=None,  # Don't save individual panels
                linewidth=[LINE_WIDTH, LINE_WIDTH],
                color=['#D62728', '#1F77B4']
            )
            return None, None  # No colorbar needed
            
        else:
            # Fat band - use temporary figure to avoid colorbar issues
            temp_fig, temp_ax = plt.subplots(figsize=(8, 6))
            
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
                vmax=vmax,
                ax=temp_ax,
                show=False,
                savefig=None,
                plot_color_bar=False,  # No colorbar in temp figure
                linewidth=[LINE_WIDTH, LINE_WIDTH]
            )
            
            # Copy all artists from temp_ax to target ax
            for line in temp_ax.get_lines():
                ax.add_line(plt.Line2D(line.get_xdata(), line.get_ydata(),
                                       color=line.get_color(),
                                       linewidth=line.get_linewidth(),
                                       linestyle=line.get_linestyle(),
                                       alpha=line.get_alpha()))
            
            # Copy collections (scatter points for fat bands)
            for collection in temp_ax.collections:
                ax.add_collection(collection)
                collection.remove()  # Remove from temp_ax
            
            plt.close(temp_fig)
            
            return cmap, vmax
            
    except Exception as e:
        print(f"    ⚠ Error: {e}")
        return None, None

def create_robust_multipanel(fermi=None):
    """
    Create multi-panel figure using robust separate plotting
    """
    setup_matplotlib()
    
    if fermi is None:
        fermi = extract_fermi_from_outcar()
        if fermi is None:
            print("\n✗ Error: Cannot determine Fermi energy")
            sys.exit(1)
    
    print(f"\nUsing Fermi energy: {fermi:.4f} eV")
    print(f"Energy window: [{ENERGY_MIN}, {ENERGY_MAX}] eV")
    
    print("\n" + "="*70)
    print("Creating robust multi-panel fat band structure")
    print("="*70)
    
    # Create figure
    fig = plt.figure(figsize=(FIGURE_WIDTH, FIGURE_HEIGHT))
    gs = GridSpec(1, 5, figure=fig, wspace=0.08, 
                  left=0.05, right=0.98, top=0.92, bottom=0.18)
    
    axes = [fig.add_subplot(gs[0, i]) for i in range(5)]
    
    # Panel configurations
    panels = [
        ('Total', 'Band Structure', 'plain', None, None, 'viridis', 1.0),
        ('Mo d orbitals', 'MoSSe Layer', 'parametric', ATOMS_MO, ORBITALS_D, COLORS['Mo'], VMAX_METAL),
        ('W d orbitals', 'WSSe Layer', 'parametric', ATOMS_W, ORBITALS_D, COLORS['W'], VMAX_METAL),
        ('S p orbitals', 'Chalcogen', 'parametric', ATOMS_S, ORBITALS_P, COLORS['S'], VMAX_CHALCOGEN),
        ('Se p orbitals', 'Chalcogen', 'parametric', ATOMS_SE, ORBITALS_P, COLORS['Se'], VMAX_CHALCOGEN),
    ]
    
    # Plot each panel separately
    for i, (title, subtitle, mode, atoms, orbitals, cmap, vmax) in enumerate(panels):
        print(f"\n  Panel {i+1}/5: {title}")
        
        ax = axes[i]
        
        # Plot the panel
        cmap_result, vmax_result = plot_single_panel(
            ax, title, subtitle, mode, atoms, orbitals, cmap, vmax, fermi, 
            show_ylabel=(i==0)
        )
        
        # Set limits and styling
        ax.set_xlim([0, KPATH_TICKS[-1]])
        ax.set_ylim([ENERGY_MIN, ENERGY_MAX])
        
        # Title
        ax.set_title(f"{title}\n{subtitle}", fontsize=TITLE_FONT_SIZE,
                    fontweight='bold', pad=12, linespacing=1.3)
        
        # Fermi level
        ax.axhline(y=0, color='black', linestyle='--', 
                  linewidth=1.5, alpha=0.8, zorder=10)
        
        # Axes
        if i == 0:
            ax.set_ylabel('$E - E_\\mathrm{F}$ (eV)', 
                         fontsize=FONT_SIZE+2, fontweight='bold')
            yticks = np.arange(ENERGY_MIN, ENERGY_MAX + 0.5, 1.0)
            ax.set_yticks(yticks)
        else:
            ax.set_yticklabels([])
        
        if i == 2:
            ax.set_xlabel('Wave Vector', fontsize=FONT_SIZE+2, fontweight='bold')
        
        # K-point labels
        ax.set_xticks(KPATH_TICKS)
        if i == 0:
            ax.set_xticklabels(KPATH_LABELS)
        else:
            ax.set_xticklabels([])
        
        # Grid
        ax.grid(True, alpha=0.25, linestyle=':', linewidth=0.8, color='gray', zorder=0)
        
        # Add colorbar for parametric modes
        if cmap_result is not None and vmax_result is not None:
            cbar_ax = fig.add_axes([
                ax.get_position().x0 + 0.005,
                0.08,
                ax.get_position().width - 0.01,
                0.025
            ])
            
            norm = mpl.colors.Normalize(vmin=0, vmax=vmax_result)
            sm = plt.cm.ScalarMappable(cmap=cmap_result, norm=norm)
            sm.set_array([])
            
            cbar = plt.colorbar(sm, cax=cbar_ax, orientation='horizontal')
            
            if vmax_result == VMAX_METAL:
                cbar.set_label('Projection Weight', fontsize=FONT_SIZE-1, fontweight='bold')
                cbar.set_ticks([0, 0.25, 0.5, 0.75, 1.0])
            else:
                cbar.set_label('Weight (×3 amplified)', fontsize=FONT_SIZE-1, fontweight='bold')
                cbar.set_ticks([0, 0.1, 0.2, 0.3])
            
            cbar.ax.tick_params(labelsize=FONT_SIZE-2)
            cbar.outline.set_linewidth(1.5)
        
        print(f"    ✓ Panel {i+1} complete")
    
    # Overall title
    fig.suptitle('MoSSe/WSSe Heterostructure: Orbital-Resolved Band Structure',
                 fontsize=TITLE_FONT_SIZE+2, fontweight='bold', y=0.98)
    
    # Save
    outputs = [
        'fatband_multipanel_robust.png',
        'fatband_multipanel_robust.pdf',
        'fatband_multipanel_robust.svg'
    ]
    
    print("\n" + "="*70)
    print("Saving figures...")
    print("="*70)
    
    for output in outputs:
        try:
            plt.savefig(output, dpi=DPI, bbox_inches='tight',
                       facecolor='white', edgecolor='none')
            size_mb = os.path.getsize(output) / (1024*1024)
            print(f"  ✓ Saved: {output} ({size_mb:.2f} MB)")
        except Exception as e:
            print(f"  ✗ Failed: {e}")
    
    plt.close()
    print("\n✓ Robust multi-panel generation complete!")

# ============================================================================
# SIMPLE 3-PANEL VERSION (Most reliable)
# ============================================================================

def create_simple_3panel(fermi=None):
    """
    Simplified 3-panel version - most reliable
    """
    if fermi is None:
        fermi = extract_fermi_from_outcar()
        if fermi is None:
            fermi = FERMI_ENERGY
    
    print("\n" + "="*70)
    print("Creating simple 3-panel version (Total + Mo d + W d)")
    print("="*70)
    
    setup_matplotlib()
    
    # Create 3 separate figures and combine
    fig, axes = plt.subplots(1, 3, figsize=(18, 7), sharey=True)
    
    configs = [
        ('Total', 'plain', None, None),
        ('Mo d', 'parametric', ATOMS_MO, ORBITALS_D),
        ('W d', 'parametric', ATOMS_W, ORBITALS_D),
    ]
    
    for i, (title, mode, atoms, orbitals) in enumerate(configs):
        print(f"\n  Panel {i+1}/3: {title}")
        
        try:
            if mode == 'plain':
                pyprocar.bandsplot(
                    code='vasp',
                    mode='plain',
                    fermi=fermi,
                    elimit=[ENERGY_MIN, ENERGY_MAX],
                    kticks=KPATH_TICKS,
                    knames=KPATH_LABELS,
                    ax=axes[i],
                    show=False,
                    savefig=None,
                    linewidth=[LINE_WIDTH, LINE_WIDTH],
                    color=['red', 'blue']
                )
            else:
                # Use savefig approach to avoid ax issues
                temp_name = f'_temp_panel_{i}.png'
                pyprocar.bandsplot(
                    code='vasp',
                    mode='parametric',
                    atoms=atoms,
                    orbitals=orbitals,
                    fermi=fermi,
                    elimit=[ENERGY_MIN, ENERGY_MAX],
                    kticks=KPATH_TICKS,
                    knames=KPATH_LABELS,
                    cmap=COLORS['Mo'] if i==1 else COLORS['W'],
                    vmin=0,
                    vmax=1,
                    savefig=temp_name,
                    show=False,
                    linewidth=[LINE_WIDTH, LINE_WIDTH]
                )
                
                # Load and display in subplot
                if os.path.exists(temp_name):
                    img = plt.imread(temp_name)
                    axes[i].imshow(img)
                    axes[i].axis('off')
                    os.remove(temp_name)
            
            axes[i].set_title(title, fontsize=18, fontweight='bold')
            print(f"    ✓ Complete")
            
        except Exception as e:
            print(f"    ✗ Error: {e}")
    
    plt.tight_layout()
    plt.savefig('fatband_3panel_simple.png', dpi=DPI, bbox_inches='tight')
    plt.savefig('fatband_3panel_simple.pdf', bbox_inches='tight')
    plt.close()
    
    print("\n✓ Simple 3-panel complete!")

# ============================================================================
# MAIN
# ============================================================================

if __name__ == '__main__':
    print("\n" + "#"*70)
    print("# Robust Multi-Panel Fat Band Structure Generator")
    print("# (Bug-resistant version)")
    print("#"*70)
    
    if not check_files():
        sys.exit(1)
    
    # Try the robust version
    try:
        create_robust_multipanel(fermi=FERMI_ENERGY)
    except Exception as e:
        print(f"\n⚠ Robust version failed: {e}")
        print("Falling back to simple version...")
        create_simple_3panel(fermi=FERMI_ENERGY)
    
    print("\n" + "#"*70)
    print("# Done!")
    print("#"*70)
