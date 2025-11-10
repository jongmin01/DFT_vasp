import os
import numpy as np
from ase.io import read
from pymatgen.io.vasp import Vasprun

# --- Configuration (Setting) ---
VRUN_PATH = "vasprun.xml"
if not os.path.exists(VRUN_PATH):
    raise FileNotFoundError(f"File {VRUN_PATH} not found in path. Please check.")

# --- Load final structure with ASE ---
# Read the last geometry (Relaxation result) from vasprun.xml.
print(f"Loading {VRUN_PATH} to extract final structure.")
try:
    atoms = read(VRUN_PATH, index=-1) # Load the last structure
except Exception as e:
    print(f"ASE loading error: {e}. Please check the VASP file again.")
    exit()

# --- 1. Lattice Constant and Cell Info ---
cell = atoms.get_cell()
lattice_constants = atoms.get_cell_lengths_and_angles()

a, b, c = lattice_constants[:3]
alpha, beta, gamma = lattice_constants[3:]

print("\n--- 1. Optimized Lattice Constants ---")
print(f"Lattice Vector A (a): {a:.4f} Å")
print(f"Lattice Vector B (b): {b:.4f} Å")
print(f"Lattice Vector C (c): {c:.4f} Å (Z-axis direction)")
print(f"Angles (α, β, γ): {alpha:.2f}°, {beta:.2f}°, {gamma:.2f}°")
# For MoS2 monolayer, A=B, alpha=beta=90, gamma=120 are ideal.

# --- 2. Calculate Mo-S Bond Length ---
mo_indices = [atom.index for atom in atoms if atom.symbol == 'Mo']
s_indices = [atom.index for atom in atoms if atom.symbol == 'S']

mo_s_bond_lengths = []
if mo_indices and s_indices:
    # Calculate the distance between one Mo atom and its 6 nearest S atoms.
    # (Since it is a monolayer structure, one Mo atom bonds with 6 S atoms, 3 above and 3 below)
    mo_pos = atoms.get_positions()[mo_indices[0]] # Position of the first Mo atom
    
    # Calculate distance considering periodic boundary conditions
    for s_index in s_indices:
        s_pos = atoms.get_positions()[s_index]
        # Set mic=True to calculate the minimum distance considering periodic boundary conditions
        dist = atoms.get_distance(mo_indices[0], s_index, mic=True)
        mo_s_bond_lengths.append(dist)

    # Remove duplicates and calculate average (bond lengths should all be identical)
    unique_bonds = np.unique(np.round(mo_s_bond_lengths, decimals=4))
    
    print("\n--- 2. Mo-S Bond Length ---")
    if len(unique_bonds) == 1:
        print(f"Average Mo-S Bond Length: {unique_bonds[0]:.4f} Å")
    else:
        # In case of various bond lengths (e.g., distorted structure)
        print(f"Unique Mo-S Bond Lengths found: {unique_bonds} Å")
        print(f"Average Mo-S Bond Length: {np.mean(mo_s_bond_lengths):.4f} Å")

else:
    print("\n--- WARNING: Could not find Mo or S atoms. Please check the atomic symbols. ---")
