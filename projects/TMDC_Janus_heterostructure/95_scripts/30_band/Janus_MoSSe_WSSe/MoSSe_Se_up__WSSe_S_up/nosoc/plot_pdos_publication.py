#!/usr/bin/env python3
"""
Publication-Quality PDOS Plotting using PyProcar
Optimized for academic papers with clean, professional appearance
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
ENERGY_MIN = -4.0  # 4 eV below Fermi
ENERGY_MAX = 4.0   # 4 eV above Fermi

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
LINEWIDTH = 1.5  # Line width

# ============================================================================
# MAIN PLOTTING FUNCTIONS
# ============================================================================

def plot_total_dos():
    """
    Plot total density of states
    """
    print("\n" + "="*60)
    print("1. Plotting Total DOS")
    print("="*60)

    try:
        pyprocar.dosplot(
            code='vasp',
            dirname='.',
            mode='plain',
            fermi=FERMI_ENERGY,
            elimit=[ENERGY_MIN, ENERGY_MAX],
            savefig='publication_dos_total.png',
            title='',
            figure_size=FIGURE_SIZE,
            dpi=DPI
        )
        print("✓ Saved: publication_dos_total.png")

        # Vector format
        pyprocar.dosplot(
            code='vasp',
            dirname='.',
            mode='plain',
            fermi=FERMI_ENERGY,
            elimit=[ENERGY_MIN, ENERGY_MAX],
            savefig='publication_dos_total.svg',
            title='',
            figure_size=FIGURE_SIZE
        )
        print("✓ Saved: publication_dos_total.svg")
    except Exception as e:
        print(f"✗ Failed: {e}")

def plot_element_resolved_dos():
    """
    Plot element-resolved PDOS
    """
    print("\n" + "="*60)
    print("2. Plotting Element-Resolved PDOS")
    print("="*60)

    elements = {
        'Mo': ATOMS_MO,
        'W': ATOMS_W,
        'S': ATOMS_S,
        'Se': ATOMS_SE
    }

    # Individual element DOS
    for element, atoms in elements.items():
        print(f"\n  → Plotting {element} PDOS...")
        try:
            pyprocar.dosplot(
                code='vasp',
                dirname='.',
                mode='parametric',
                atoms=atoms,
                fermi=FERMI_ENERGY,
                elimit=[ENERGY_MIN, ENERGY_MAX],
                savefig=f'publication_pdos_{element}.png',
                title='',
                figure_size=FIGURE_SIZE,
                dpi=DPI,
                labels=[element]
            )
            print(f"    ✓ Saved: publication_pdos_{element}.png")
        except Exception as e:
            print(f"    ✗ Failed: {e}")

    # Combined element DOS
    print(f"\n  → Plotting combined element PDOS...")
    try:
        pyprocar.dosplot(
            code='vasp',
            dirname='.',
            mode='stack',
            atoms=[ATOMS_MO, ATOMS_W, ATOMS_S, ATOMS_SE],
            fermi=FERMI_ENERGY,
            elimit=[ENERGY_MIN, ENERGY_MAX],
            savefig='publication_pdos_elements_stacked.png',
            title='',
            figure_size=FIGURE_SIZE,
            dpi=DPI,
            labels=['Mo', 'W', 'S', 'Se']
        )
        print(f"    ✓ Saved: publication_pdos_elements_stacked.png")
    except Exception as e:
        print(f"    Note: Stacked mode not available, trying overlay...")
        try:
            pyprocar.dosplot(
                code='vasp',
                dirname='.',
                mode='parametric_line',
                atoms=[ATOMS_MO, ATOMS_W, ATOMS_S, ATOMS_SE],
                fermi=FERMI_ENERGY,
                elimit=[ENERGY_MIN, ENERGY_MAX],
                savefig='publication_pdos_elements_overlay.png',
                title='',
                figure_size=FIGURE_SIZE,
                dpi=DPI,
                labels=['Mo', 'W', 'S', 'Se']
            )
            print(f"    ✓ Saved: publication_pdos_elements_overlay.png")
        except Exception as e2:
            print(f"    ✗ Failed: {e2}")

def plot_layer_resolved_dos():
    """
    Plot layer-resolved PDOS (MoSSe vs WSSe)
    """
    print("\n" + "="*60)
    print("3. Plotting Layer-Resolved PDOS")
    print("="*60)

    # MoSSe layer
    print("\n  → Plotting MoSSe layer PDOS...")
    try:
        pyprocar.dosplot(
            code='vasp',
            dirname='.',
            mode='parametric',
            atoms=ATOMS_MO,
            fermi=FERMI_ENERGY,
            elimit=[ENERGY_MIN, ENERGY_MAX],
            savefig='publication_pdos_MoSSe_layer.png',
            title='',
            figure_size=FIGURE_SIZE,
            dpi=DPI,
            labels=['MoSSe']
        )
        print("    ✓ Saved: publication_pdos_MoSSe_layer.png")
    except Exception as e:
        print(f"    ✗ Failed: {e}")

    # WSSe layer
    print("\n  → Plotting WSSe layer PDOS...")
    try:
        pyprocar.dosplot(
            code='vasp',
            dirname='.',
            mode='parametric',
            atoms=ATOMS_W,
            fermi=FERMI_ENERGY,
            elimit=[ENERGY_MIN, ENERGY_MAX],
            savefig='publication_pdos_WSSe_layer.png',
            title='',
            figure_size=FIGURE_SIZE,
            dpi=DPI,
            labels=['WSSe']
        )
        print("    ✓ Saved: publication_pdos_WSSe_layer.png")
    except Exception as e:
        print(f"    ✗ Failed: {e}")

    # Layer comparison
    print("\n  → Plotting layer comparison...")
    try:
        pyprocar.dosplot(
            code='vasp',
            dirname='.',
            mode='parametric_line',
            atoms=[ATOMS_MO, ATOMS_W],
            fermi=FERMI_ENERGY,
            elimit=[ENERGY_MIN, ENERGY_MAX],
            savefig='publication_pdos_layer_comparison.png',
            title='',
            figure_size=FIGURE_SIZE,
            dpi=DPI,
            labels=['MoSSe', 'WSSe']
        )
        print("    ✓ Saved: publication_pdos_layer_comparison.png")
    except Exception as e:
        print(f"    ✗ Failed: {e}")

def plot_orbital_resolved_dos():
    """
    Plot orbital-resolved PDOS
    """
    print("\n" + "="*60)
    print("4. Plotting Orbital-Resolved PDOS")
    print("="*60)

    # Mo d-orbitals
    print("\n  → Plotting Mo d-orbital PDOS...")
    try:
        pyprocar.dosplot(
            code='vasp',
            dirname='.',
            mode='parametric',
            atoms=ATOMS_MO,
            orbitals=[4, 5, 6, 7, 8],  # d-orbitals
            fermi=FERMI_ENERGY,
            elimit=[ENERGY_MIN, ENERGY_MAX],
            savefig='publication_pdos_Mo_d.png',
            title='',
            figure_size=FIGURE_SIZE,
            dpi=DPI,
            labels=['Mo d']
        )
        print("    ✓ Saved: publication_pdos_Mo_d.png")
    except Exception as e:
        print(f"    ✗ Failed: {e}")

    # W d-orbitals
    print("\n  → Plotting W d-orbital PDOS...")
    try:
        pyprocar.dosplot(
            code='vasp',
            dirname='.',
            mode='parametric',
            atoms=ATOMS_W,
            orbitals=[4, 5, 6, 7, 8],  # d-orbitals
            fermi=FERMI_ENERGY,
            elimit=[ENERGY_MIN, ENERGY_MAX],
            savefig='publication_pdos_W_d.png',
            title='',
            figure_size=FIGURE_SIZE,
            dpi=DPI,
            labels=['W d']
        )
        print("    ✓ Saved: publication_pdos_W_d.png")
    except Exception as e:
        print(f"    ✗ Failed: {e}")

    # S p-orbitals
    print("\n  → Plotting S p-orbital PDOS...")
    try:
        pyprocar.dosplot(
            code='vasp',
            dirname='.',
            mode='parametric',
            atoms=ATOMS_S,
            orbitals=[1, 2, 3],  # p-orbitals
            fermi=FERMI_ENERGY,
            elimit=[ENERGY_MIN, ENERGY_MAX],
            savefig='publication_pdos_S_p.png',
            title='',
            figure_size=FIGURE_SIZE,
            dpi=DPI,
            labels=['S p']
        )
        print("    ✓ Saved: publication_pdos_S_p.png")
    except Exception as e:
        print(f"    ✗ Failed: {e}")

    # Se p-orbitals
    print("\n  → Plotting Se p-orbital PDOS...")
    try:
        pyprocar.dosplot(
            code='vasp',
            dirname='.',
            mode='parametric',
            atoms=ATOMS_SE,
            orbitals=[1, 2, 3],  # p-orbitals
            fermi=FERMI_ENERGY,
            elimit=[ENERGY_MIN, ENERGY_MAX],
            savefig='publication_pdos_Se_p.png',
            title='',
            figure_size=FIGURE_SIZE,
            dpi=DPI,
            labels=['Se p']
        )
        print("    ✓ Saved: publication_pdos_Se_p.png")
    except Exception as e:
        print(f"    ✗ Failed: {e}")

def plot_chalcogen_comparison():
    """
    Plot S vs Se comparison
    """
    print("\n" + "="*60)
    print("5. Plotting Chalcogen Comparison (S vs Se)")
    print("="*60)

    try:
        pyprocar.dosplot(
            code='vasp',
            dirname='.',
            mode='parametric_line',
            atoms=[ATOMS_S, ATOMS_SE],
            fermi=FERMI_ENERGY,
            elimit=[ENERGY_MIN, ENERGY_MAX],
            savefig='publication_pdos_S_vs_Se.png',
            title='',
            figure_size=FIGURE_SIZE,
            dpi=DPI,
            labels=['S', 'Se']
        )
        print("✓ Saved: publication_pdos_S_vs_Se.png")
    except Exception as e:
        print(f"✗ Failed: {e}")

def plot_metal_comparison():
    """
    Plot Mo vs W comparison
    """
    print("\n" + "="*60)
    print("6. Plotting Metal Comparison (Mo vs W)")
    print("="*60)

    try:
        pyprocar.dosplot(
            code='vasp',
            dirname='.',
            mode='parametric_line',
            atoms=[ATOMS_MO, ATOMS_W],
            fermi=FERMI_ENERGY,
            elimit=[ENERGY_MIN, ENERGY_MAX],
            savefig='publication_pdos_Mo_vs_W.png',
            title='',
            figure_size=FIGURE_SIZE,
            dpi=DPI,
            labels=['Mo', 'W']
        )
        print("✓ Saved: publication_pdos_Mo_vs_W.png")
    except Exception as e:
        print(f"✗ Failed: {e}")

# ============================================================================
# MAIN EXECUTION
# ============================================================================

if __name__ == '__main__':
    print("\n" + "#"*60)
    print("# Publication-Quality PDOS Generation")
    print("#"*60)
    print(f"\nConfiguration:")
    print(f"  Fermi Energy: {FERMI_ENERGY} eV")
    print(f"  Energy Range: [{ENERGY_MIN}, {ENERGY_MAX}] eV (relative to Ef)")
    print(f"  Resolution: {DPI} DPI")
    print(f"  Figure Size: {FIGURE_SIZE} inches")

    # Generate plots
    plot_total_dos()
    plot_element_resolved_dos()
    plot_layer_resolved_dos()
    plot_orbital_resolved_dos()
    plot_chalcogen_comparison()
    plot_metal_comparison()

    print("\n" + "#"*60)
    print("# All Publication-Quality PDOS Plots Generated!")
    print("#"*60)
    print("\nGenerated files:")
    print("  Total DOS:")
    print("    - publication_dos_total.png/svg")
    print("\n  Element-resolved:")
    print("    - publication_pdos_Mo.png")
    print("    - publication_pdos_W.png")
    print("    - publication_pdos_S.png")
    print("    - publication_pdos_Se.png")
    print("    - publication_pdos_elements_*.png")
    print("\n  Layer-resolved:")
    print("    - publication_pdos_MoSSe_layer.png")
    print("    - publication_pdos_WSSe_layer.png")
    print("    - publication_pdos_layer_comparison.png")
    print("\n  Orbital-resolved:")
    print("    - publication_pdos_Mo_d.png")
    print("    - publication_pdos_W_d.png")
    print("    - publication_pdos_S_p.png")
    print("    - publication_pdos_Se_p.png")
    print("\n  Comparisons:")
    print("    - publication_pdos_S_vs_Se.png")
    print("    - publication_pdos_Mo_vs_W.png")
    print("\n" + "#"*60)
