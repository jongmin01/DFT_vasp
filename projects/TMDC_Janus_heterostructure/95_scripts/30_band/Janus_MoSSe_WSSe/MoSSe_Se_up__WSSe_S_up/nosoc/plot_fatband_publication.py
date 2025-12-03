#!/usr/bin/env python3
"""
Publication-Quality Fat Band Structure Plotting using PyProcar
Optimized for academic papers with Fermi-level centered energy range
"""

import pyprocar
import matplotlib.pyplot as plt
import numpy as np

# ============================================================================
# CONFIGURATION - Adjust these for your needs
# ============================================================================

# Fermi energy (eV) - obtained from OUTCAR
FERMI_ENERGY = -0.302

# Energy range relative to Fermi level (eV)
ENERGY_MIN = -3.0  # 3 eV below Fermi
ENERGY_MAX = 3.0   # 3 eV above Fermi

# High-symmetry k-points for Γ-M-K-Γ path
KPATH_LABELS = ['$\\Gamma$', 'M', 'K', '$\\Gamma$']
KPATH_TICKS = [0, 43, 83, 120]

# Atom indices (0-indexed)
# From POSCAR: Mo(0-24), W(25-49), S(50-99), Se(100-149)
ATOMS_MO = list(range(0, 25))
ATOMS_W = list(range(25, 50))
ATOMS_S = list(range(50, 100))
ATOMS_SE = list(range(100, 150))

# Publication settings
DPI = 600  # High resolution for publication
FIGURE_SIZE = (8, 6)  # Figure size in inches
FONT_SIZE = 14  # Base font size
LINEWIDTH = (1.5, 1.5)  # Line width for bands (spin up, spin down)

# Color schemes for different elements
COLORS = {
    'Mo': 'Reds',
    'W': 'Blues',
    'S': 'Greens',
    'Se': 'Purples'
}

# ============================================================================
# MAIN PLOTTING FUNCTIONS
# ============================================================================

def plot_total_band_structure():
    """
    Plot clean total band structure without projections
    """
    print("\n" + "="*60)
    print("1. Plotting Total Band Structure")
    print("="*60)

    pyprocar.bandsplot(
        code='vasp',
        dirname='.',
        mode='plain',
        fermi=FERMI_ENERGY,
        elimit=[ENERGY_MIN, ENERGY_MAX],
        kticks=KPATH_TICKS,
        knames=KPATH_LABELS,
        savefig='publication_band_total.png',
        title='',  # No title for publication
        figure_size=FIGURE_SIZE,
        dpi=DPI,
        linewidth=LINEWIDTH
    )
    print("✓ Saved: publication_band_total.png")

    # Also save as vector format
    try:
        pyprocar.bandsplot(
            code='vasp',
            dirname='.',
            mode='plain',
            fermi=FERMI_ENERGY,
            elimit=[ENERGY_MIN, ENERGY_MAX],
            kticks=KPATH_TICKS,
            knames=KPATH_LABELS,
            savefig='publication_band_total.svg',
            title='',
            figure_size=FIGURE_SIZE,
            linewidth=LINEWIDTH
        )
        print("✓ Saved: publication_band_total.svg (vector format)")
    except:
        pass

def plot_element_projected_bands():
    """
    Plot element-resolved fat band structures
    """
    print("\n" + "="*60)
    print("2. Plotting Element-Projected Band Structures")
    print("="*60)

    elements = {
        'Mo': (ATOMS_MO, COLORS['Mo']),
        'W': (ATOMS_W, COLORS['W']),
        'S': (ATOMS_S, COLORS['S']),
        'Se': (ATOMS_SE, COLORS['Se'])
    }

    for element, (atoms, cmap) in elements.items():
        print(f"\n  → Plotting {element} contribution...")

        pyprocar.bandsplot(
            code='vasp',
            dirname='.',
            mode='parametric',
            atoms=atoms,
            fermi=FERMI_ENERGY,
            elimit=[ENERGY_MIN, ENERGY_MAX],
            kticks=KPATH_TICKS,
            knames=KPATH_LABELS,
            cmap=cmap,
            vmin=0,
            vmax=1,
            savefig=f'publication_band_{element}.png',
            title='',
            figure_size=FIGURE_SIZE,
            dpi=DPI,
            linewidth=LINEWIDTH
        )
        print(f"    ✓ Saved: publication_band_{element}.png")

def plot_layer_comparison():
    """
    Plot MoSSe vs WSSe layer comparison
    """
    print("\n" + "="*60)
    print("3. Plotting Layer-Resolved Band Structures")
    print("="*60)

    # MoSSe layer (Mo atoms)
    print("\n  → Plotting MoSSe layer...")
    pyprocar.bandsplot(
        code='vasp',
        dirname='.',
        mode='parametric',
        atoms=ATOMS_MO,
        fermi=FERMI_ENERGY,
        elimit=[ENERGY_MIN, ENERGY_MAX],
        kticks=KPATH_TICKS,
        knames=KPATH_LABELS,
        cmap='Reds',
        vmin=0,
        vmax=1,
        savefig='publication_band_MoSSe_layer.png',
        title='',
        figure_size=FIGURE_SIZE,
        dpi=DPI,
        linewidth=LINEWIDTH
    )
    print("    ✓ Saved: publication_band_MoSSe_layer.png")

    # WSSe layer (W atoms)
    print("\n  → Plotting WSSe layer...")
    pyprocar.bandsplot(
        code='vasp',
        dirname='.',
        mode='parametric',
        atoms=ATOMS_W,
        fermi=FERMI_ENERGY,
        elimit=[ENERGY_MIN, ENERGY_MAX],
        kticks=KPATH_TICKS,
        knames=KPATH_LABELS,
        cmap='Blues',
        vmin=0,
        vmax=1,
        savefig='publication_band_WSSe_layer.png',
        title='',
        figure_size=FIGURE_SIZE,
        dpi=DPI,
        linewidth=LINEWIDTH
    )
    print("    ✓ Saved: publication_band_WSSe_layer.png")

def plot_orbital_resolved():
    """
    Plot orbital-resolved band structures for transition metals
    """
    print("\n" + "="*60)
    print("4. Plotting Orbital-Resolved Band Structures")
    print("="*60)

    # d-orbitals for Mo
    print("\n  → Plotting Mo d-orbitals...")
    try:
        pyprocar.bandsplot(
            code='vasp',
            dirname='.',
            mode='parametric',
            atoms=ATOMS_MO,
            orbitals=[4, 5, 6, 7, 8],  # d-orbitals
            fermi=FERMI_ENERGY,
            elimit=[ENERGY_MIN, ENERGY_MAX],
            kticks=KPATH_TICKS,
            knames=KPATH_LABELS,
            cmap='Reds',
            vmin=0,
            vmax=1,
            savefig='publication_band_Mo_d.png',
            title='',
            figure_size=FIGURE_SIZE,
            dpi=DPI,
            linewidth=LINEWIDTH
        )
        print("    ✓ Saved: publication_band_Mo_d.png")
    except Exception as e:
        print(f"    ✗ Failed: {e}")

    # d-orbitals for W
    print("\n  → Plotting W d-orbitals...")
    try:
        pyprocar.bandsplot(
            code='vasp',
            dirname='.',
            mode='parametric',
            atoms=ATOMS_W,
            orbitals=[4, 5, 6, 7, 8],  # d-orbitals
            fermi=FERMI_ENERGY,
            elimit=[ENERGY_MIN, ENERGY_MAX],
            kticks=KPATH_TICKS,
            knames=KPATH_LABELS,
            cmap='Blues',
            vmin=0,
            vmax=1,
            savefig='publication_band_W_d.png',
            title='',
            figure_size=FIGURE_SIZE,
            dpi=DPI,
            linewidth=LINEWIDTH
        )
        print("    ✓ Saved: publication_band_W_d.png")
    except Exception as e:
        print(f"    ✗ Failed: {e}")

# ============================================================================
# MAIN EXECUTION
# ============================================================================

if __name__ == '__main__':
    print("\n" + "#"*60)
    print("# Publication-Quality Fat Band Structure Generation")
    print("#"*60)
    print(f"\nConfiguration:")
    print(f"  Fermi Energy: {FERMI_ENERGY} eV")
    print(f"  Energy Range: [{ENERGY_MIN}, {ENERGY_MAX}] eV (relative to Ef)")
    print(f"  Resolution: {DPI} DPI")
    print(f"  Figure Size: {FIGURE_SIZE} inches")
    print(f"  K-path: {' → '.join(KPATH_LABELS)}")

    # Generate plots
    plot_total_band_structure()
    plot_element_projected_bands()
    plot_layer_comparison()
    plot_orbital_resolved()

    print("\n" + "#"*60)
    print("# All Publication-Quality Band Structure Plots Generated!")
    print("#"*60)
    print("\nGenerated files:")
    print("  Main plots:")
    print("    - publication_band_total.png/svg")
    print("    - publication_band_Mo.png")
    print("    - publication_band_W.png")
    print("    - publication_band_S.png")
    print("    - publication_band_Se.png")
    print("\n  Layer-resolved:")
    print("    - publication_band_MoSSe_layer.png")
    print("    - publication_band_WSSe_layer.png")
    print("\n  Orbital-resolved:")
    print("    - publication_band_Mo_d.png")
    print("    - publication_band_W_d.png")
    print("\n" + "#"*60)
