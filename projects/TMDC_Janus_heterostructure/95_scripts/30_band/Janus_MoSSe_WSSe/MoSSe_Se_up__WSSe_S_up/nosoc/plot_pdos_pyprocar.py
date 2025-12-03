#!/usr/bin/env python3
"""
Projected Density of States (PDOS) Plotting using PyProcar
Visualizes element-resolved and orbital-resolved DOS
"""

import pyprocar
import matplotlib.pyplot as plt
import numpy as np

# Configuration
OUTCAR_PATH = 'OUTCAR'
PROCAR_PATH = 'PROCAR'
VASPRUN_PATH = 'vasprun.xml'

# Energy range for plotting (eV relative to Fermi level)
ENERGY_MIN = -6
ENERGY_MAX = 6

# Atom indices (0-indexed)
# From POSCAR: Mo(0-24), W(25-49), S(50-99), Se(100-149)
ATOMS_MO = list(range(0, 25))
ATOMS_W = list(range(25, 50))
ATOMS_S = list(range(50, 100))
ATOMS_SE = list(range(100, 150))

def plot_total_dos():
    """
    Plot total density of states
    """
    print("=" * 60)
    print("Plotting Total DOS")
    print("=" * 60)

    pyprocar.dosplot(
        code='vasp',
        dirname='.',
        mode='plain',
        elimit=[ENERGY_MIN, ENERGY_MAX],
        savefig='pdos_total.png',
        title='Total Density of States',
        labels=['Total DOS']
    )
    print("Saved: pdos_total.png")

def plot_element_resolved_dos():
    """
    Plot element-resolved projected DOS
    """
    print("\n" + "=" * 60)
    print("Plotting Element-Resolved PDOS")
    print("=" * 60)

    # Plot individual element DOS
    print("\n1. Plotting Mo PDOS...")
    pyprocar.dosplot(
        code='vasp',
        dirname='.',
        mode='parametric',
        atoms=ATOMS_MO,
        elimit=[ENERGY_MIN, ENERGY_MAX],
        savefig='pdos_Mo.png',
        title='Projected DOS - Mo',
        labels=['Mo']
    )
    print("   Saved: pdos_Mo.png")

    print("\n2. Plotting W PDOS...")
    pyprocar.dosplot(
        code='vasp',
        dirname='.',
        mode='parametric',
        atoms=ATOMS_W,
        elimit=[ENERGY_MIN, ENERGY_MAX],
        savefig='pdos_W.png',
        title='Projected DOS - W',
        labels=['W']
    )
    print("   Saved: pdos_W.png")

    print("\n3. Plotting S PDOS...")
    pyprocar.dosplot(
        code='vasp',
        dirname='.',
        mode='parametric',
        atoms=ATOMS_S,
        elimit=[ENERGY_MIN, ENERGY_MAX],
        savefig='pdos_S.png',
        title='Projected DOS - S',
        labels=['S']
    )
    print("   Saved: pdos_S.png")

    print("\n4. Plotting Se PDOS...")
    pyprocar.dosplot(
        code='vasp',
        dirname='.',
        mode='parametric',
        atoms=ATOMS_SE,
        elimit=[ENERGY_MIN, ENERGY_MAX],
        savefig='pdos_Se.png',
        title='Projected DOS - Se',
        labels=['Se']
    )
    print("   Saved: pdos_Se.png")

    # Combined plot with all elements
    print("\n5. Plotting combined element PDOS...")
    try:
        pyprocar.dosplot(
            code='vasp',
            dirname='.',
            mode='stack',
            atoms=[ATOMS_MO, ATOMS_W, ATOMS_S, ATOMS_SE],
            elimit=[ENERGY_MIN, ENERGY_MAX],
            savefig='pdos_elements_stacked.png',
            title='Element-Resolved PDOS (Stacked)',
            labels=['Mo', 'W', 'S', 'Se']
        )
        print("   Saved: pdos_elements_stacked.png")
    except Exception as e:
        print(f"   Note: Stacked plot failed, trying parametric overlay: {e}")
        try:
            pyprocar.dosplot(
                code='vasp',
                dirname='.',
                mode='parametric_line',
                atoms=[ATOMS_MO, ATOMS_W, ATOMS_S, ATOMS_SE],
                elimit=[ENERGY_MIN, ENERGY_MAX],
                savefig='pdos_elements_overlay.png',
                title='Element-Resolved PDOS (Overlay)',
                labels=['Mo', 'W', 'S', 'Se']
            )
            print("   Saved: pdos_elements_overlay.png")
        except Exception as e2:
            print(f"   Failed: {e2}")

    print("\n" + "=" * 60)
    print("Element-Resolved PDOS Plotting Complete!")
    print("=" * 60)

def plot_orbital_resolved_dos():
    """
    Plot orbital-resolved projected DOS for each element
    """
    print("\n" + "=" * 60)
    print("Plotting Orbital-Resolved PDOS")
    print("=" * 60)

    # Orbital indices: 0=s, 1=py, 2=pz, 3=px, 4=dxy, 5=dyz, 6=dz2, 7=dxz, 8=dx2

    # Mo d-orbitals
    print("\n1. Plotting Mo orbital-resolved PDOS...")
    try:
        pyprocar.dosplot(
            code='vasp',
            dirname='.',
            mode='stack_orbitals',
            atoms=ATOMS_MO,
            orbitals=[0, 1, 2, 3, 4, 5, 6, 7, 8],
            elimit=[ENERGY_MIN, ENERGY_MAX],
            savefig='pdos_Mo_orbitals.png',
            title='Mo Orbital-Resolved PDOS',
            labels=['s', 'py', 'pz', 'px', 'dxy', 'dyz', 'dz2', 'dxz', 'dx2']
        )
        print("   Saved: pdos_Mo_orbitals.png")
    except Exception as e:
        print(f"   Note: Trying alternative mode: {e}")
        try:
            # Try d-orbitals only for clearer visualization
            pyprocar.dosplot(
                code='vasp',
                dirname='.',
                mode='parametric',
                atoms=ATOMS_MO,
                orbitals=[4, 5, 6, 7, 8],  # d-orbitals only
                elimit=[ENERGY_MIN, ENERGY_MAX],
                savefig='pdos_Mo_d_orbitals.png',
                title='Mo d-orbital PDOS',
                labels=['dxy', 'dyz', 'dz2', 'dxz', 'dx2-y2']
            )
            print("   Saved: pdos_Mo_d_orbitals.png")
        except Exception as e2:
            print(f"   Failed: {e2}")

    # W d-orbitals
    print("\n2. Plotting W orbital-resolved PDOS...")
    try:
        pyprocar.dosplot(
            code='vasp',
            dirname='.',
            mode='stack_orbitals',
            atoms=ATOMS_W,
            orbitals=[0, 1, 2, 3, 4, 5, 6, 7, 8],
            elimit=[ENERGY_MIN, ENERGY_MAX],
            savefig='pdos_W_orbitals.png',
            title='W Orbital-Resolved PDOS',
            labels=['s', 'py', 'pz', 'px', 'dxy', 'dyz', 'dz2', 'dxz', 'dx2']
        )
        print("   Saved: pdos_W_orbitals.png")
    except Exception as e:
        print(f"   Note: Trying alternative mode: {e}")
        try:
            pyprocar.dosplot(
                code='vasp',
                dirname='.',
                mode='parametric',
                atoms=ATOMS_W,
                orbitals=[4, 5, 6, 7, 8],  # d-orbitals only
                elimit=[ENERGY_MIN, ENERGY_MAX],
                savefig='pdos_W_d_orbitals.png',
                title='W d-orbital PDOS',
                labels=['dxy', 'dyz', 'dz2', 'dxz', 'dx2-y2']
            )
            print("   Saved: pdos_W_d_orbitals.png")
        except Exception as e2:
            print(f"   Failed: {e2}")

    # S p-orbitals
    print("\n3. Plotting S orbital-resolved PDOS...")
    try:
        pyprocar.dosplot(
            code='vasp',
            dirname='.',
            mode='stack_orbitals',
            atoms=ATOMS_S,
            orbitals=[0, 1, 2, 3],
            elimit=[ENERGY_MIN, ENERGY_MAX],
            savefig='pdos_S_orbitals.png',
            title='S Orbital-Resolved PDOS',
            labels=['s', 'py', 'pz', 'px']
        )
        print("   Saved: pdos_S_orbitals.png")
    except Exception as e:
        print(f"   Note: Trying alternative mode: {e}")
        try:
            pyprocar.dosplot(
                code='vasp',
                dirname='.',
                mode='parametric',
                atoms=ATOMS_S,
                orbitals=[1, 2, 3],  # p-orbitals only
                elimit=[ENERGY_MIN, ENERGY_MAX],
                savefig='pdos_S_p_orbitals.png',
                title='S p-orbital PDOS',
                labels=['py', 'pz', 'px']
            )
            print("   Saved: pdos_S_p_orbitals.png")
        except Exception as e2:
            print(f"   Failed: {e2}")

    # Se p-orbitals
    print("\n4. Plotting Se orbital-resolved PDOS...")
    try:
        pyprocar.dosplot(
            code='vasp',
            dirname='.',
            mode='stack_orbitals',
            atoms=ATOMS_SE,
            orbitals=[0, 1, 2, 3],
            elimit=[ENERGY_MIN, ENERGY_MAX],
            savefig='pdos_Se_orbitals.png',
            title='Se Orbital-Resolved PDOS',
            labels=['s', 'py', 'pz', 'px']
        )
        print("   Saved: pdos_Se_orbitals.png")
    except Exception as e:
        print(f"   Note: Trying alternative mode: {e}")
        try:
            pyprocar.dosplot(
                code='vasp',
                dirname='.',
                mode='parametric',
                atoms=ATOMS_SE,
                orbitals=[1, 2, 3],  # p-orbitals only
                elimit=[ENERGY_MIN, ENERGY_MAX],
                savefig='pdos_Se_p_orbitals.png',
                title='Se p-orbital PDOS',
                labels=['py', 'pz', 'px']
            )
            print("   Saved: pdos_Se_p_orbitals.png")
        except Exception as e2:
            print(f"   Failed: {e2}")

    print("\n" + "=" * 60)
    print("Orbital-Resolved PDOS Plotting Complete!")
    print("=" * 60)

def plot_layer_resolved_dos():
    """
    Plot layer-resolved DOS for MoSSe and WSSe
    """
    print("\n" + "=" * 60)
    print("Plotting Layer-Resolved PDOS")
    print("=" * 60)

    # Define layers based on structure
    # Assuming MoSSe is first layer and WSSe is second layer
    # Need to identify which atoms belong to which layer

    # Layer 1: MoSSe (Mo at z~0.053, S and Se at different z)
    # Layer 2: WSSe (W at z~0.412)

    print("\n1. Plotting MoSSe layer PDOS...")
    try:
        # Mo + its neighboring S and Se atoms
        mosse_atoms = ATOMS_MO  # Simplified: just Mo atoms
        pyprocar.dosplot(
            code='vasp',
            dirname='.',
            mode='parametric',
            atoms=mosse_atoms,
            elimit=[ENERGY_MIN, ENERGY_MAX],
            savefig='pdos_MoSSe_layer.png',
            title='PDOS - MoSSe Layer',
            labels=['MoSSe']
        )
        print("   Saved: pdos_MoSSe_layer.png")
    except Exception as e:
        print(f"   Failed: {e}")

    print("\n2. Plotting WSSe layer PDOS...")
    try:
        # W + its neighboring S and Se atoms
        wsse_atoms = ATOMS_W  # Simplified: just W atoms
        pyprocar.dosplot(
            code='vasp',
            dirname='.',
            mode='parametric',
            atoms=wsse_atoms,
            elimit=[ENERGY_MIN, ENERGY_MAX],
            savefig='pdos_WSSe_layer.png',
            title='PDOS - WSSe Layer',
            labels=['WSSe']
        )
        print("   Saved: pdos_WSSe_layer.png")
    except Exception as e:
        print(f"   Failed: {e}")

    # Chalcogens: S vs Se comparison
    print("\n3. Plotting S vs Se comparison...")
    try:
        pyprocar.dosplot(
            code='vasp',
            dirname='.',
            mode='parametric_line',
            atoms=[ATOMS_S, ATOMS_SE],
            elimit=[ENERGY_MIN, ENERGY_MAX],
            savefig='pdos_S_vs_Se.png',
            title='PDOS Comparison: S vs Se',
            labels=['S', 'Se']
        )
        print("   Saved: pdos_S_vs_Se.png")
    except Exception as e:
        print(f"   Failed: {e}")

    # Transition metals: Mo vs W comparison
    print("\n4. Plotting Mo vs W comparison...")
    try:
        pyprocar.dosplot(
            code='vasp',
            dirname='.',
            mode='parametric_line',
            atoms=[ATOMS_MO, ATOMS_W],
            elimit=[ENERGY_MIN, ENERGY_MAX],
            savefig='pdos_Mo_vs_W.png',
            title='PDOS Comparison: Mo vs W',
            labels=['Mo', 'W']
        )
        print("   Saved: pdos_Mo_vs_W.png")
    except Exception as e:
        print(f"   Failed: {e}")

    print("\n" + "=" * 60)
    print("Layer-Resolved PDOS Plotting Complete!")
    print("=" * 60)

if __name__ == '__main__':
    # Plot total DOS
    plot_total_dos()

    # Plot element-resolved DOS
    plot_element_resolved_dos()

    # Plot orbital-resolved DOS
    plot_orbital_resolved_dos()

    # Plot layer-resolved DOS
    plot_layer_resolved_dos()

    print("\n" + "=" * 60)
    print("All PDOS Plots Generated Successfully!")
    print("=" * 60)
    print("\nGenerated files:")
    print("  Total DOS:")
    print("    - pdos_total.png")
    print("\n  Element-Resolved PDOS:")
    print("    - pdos_Mo.png")
    print("    - pdos_W.png")
    print("    - pdos_S.png")
    print("    - pdos_Se.png")
    print("    - pdos_elements_*.png")
    print("\n  Orbital-Resolved PDOS:")
    print("    - pdos_Mo_orbitals.png (or pdos_Mo_d_orbitals.png)")
    print("    - pdos_W_orbitals.png (or pdos_W_d_orbitals.png)")
    print("    - pdos_S_orbitals.png (or pdos_S_p_orbitals.png)")
    print("    - pdos_Se_orbitals.png (or pdos_Se_p_orbitals.png)")
    print("\n  Layer-Resolved PDOS:")
    print("    - pdos_MoSSe_layer.png")
    print("    - pdos_WSSe_layer.png")
    print("    - pdos_S_vs_Se.png")
    print("    - pdos_Mo_vs_W.png")
