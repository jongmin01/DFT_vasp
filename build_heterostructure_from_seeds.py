#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
build_heterostructure_from_seeds.py
Build heterostructures from primitive Janus seed files
Python 2/3 compatible
"""

from __future__ import print_function
import os
import sys
import numpy as np

# ============================================================
# USER SETTINGS
# ============================================================

# Input structure files (4 Janus structures)
structures = {
    'MoSSe_S_up': './MoSeS_S_up.vasp',
    'MoSSe_Se_up': './MoSeS_Se_up.vasp',
    'WSSe_S_up': './WSeS_S_up.vasp',
    'WSSe_Se_up': './WSeS_Se_up.vasp'
}

# Structural parameters
supercell = (5, 5, 1)               # Supercell expansion
interlayer_distance = 6.8           # Interlayer distance in Angstrom (initial guess)
cell_height_c = 30.0                # Total cell height in Angstrom

# INCAR settings
incar_template = {
    'SYSTEM': 'Janus TMD Heterostructure',
    'PREC': 'Normal',
    'ENCUT': 400,
    'EDIFF': 1e-5,
    'ALGO': 'Normal',
    'NELM': 150,
    'ISMEAR': 0,
    'SIGMA': 0.1,
    'IBRION': 2,
    'ISIF': 2,
    'NSW': 200,
    'EDIFFG': -0.02,
    'AMIX': 0.2,
    'BMIX': 0.0001,
    'MAXMIX': 40,
    'LREAL': '.FALSE.',
    'ADDGRID': '.TRUE.',
    'LWAVE': '.FALSE.',
    'LCHARG': '.FALSE.',
    'NCORE': 4,
}

# All stacking combinations to calculate
combinations_to_run = [
    ('MoSSe_S_up', 'WSSe_S_up'),
    ('MoSSe_S_up', 'WSSe_Se_up'),
    ('MoSSe_Se_up', 'WSSe_S_up'),
    ('MoSSe_Se_up', 'WSSe_Se_up'),
]

# ============================================================
# FUNCTIONS
# ============================================================

def read_vasp_file(filename):
    """Read VASP POSCAR file and return structure info"""
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    comment = lines[0].strip()
    scale = float(lines[1].strip())
    
    # Lattice vectors
    lattice = []
    for i in range(2, 5):
        vec = [float(x) for x in lines[i].split()]
        lattice.append([v * scale for v in vec])
    
    # Species and counts
    species = lines[5].split()
    counts = [int(x) for x in lines[6].split()]
    
    # Coordinate type
    coord_type = lines[7].strip()
    
    # Coordinates
    coords = []
    for i in range(8, 8 + sum(counts)):
        coord = [float(x) for x in lines[i].split()[:3]]
        coords.append(coord)
    
    return {
        'comment': comment,
        'lattice': lattice,
        'species': species,
        'counts': counts,
        'coord_type': coord_type,
        'coords': coords,
    }

def make_supercell(structure, nx, ny, nz):
    """Create supercell from unit cell"""
    coords = structure['coords']
    species = structure['species']
    counts = structure['counts']
    
    super_coords = []
    super_species = []
    
    # Build atom type list
    atom_types = []
    for sp, count in zip(species, counts):
        atom_types.extend([sp] * count)
    
    # Replicate atoms
    for ix in range(nx):
        for iy in range(ny):
            for iz in range(nz):
                for coord, typ in zip(coords, atom_types):
                    new_coord = [
                        (coord[0] + ix) / float(nx),
                        (coord[1] + iy) / float(ny),
                        coord[2]  # Keep z fractional as is
                    ]
                    super_coords.append(new_coord)
                    super_species.append(typ)
    
    return super_coords, super_species

def shift_layer_z(coords, z_shift, c_lattice):
    """Shift layer in z direction"""
    shifted_coords = []
    for coord in coords:
        # Convert fractional z to Cartesian
        z_cart = coord[2] * c_lattice
        # Add shift
        z_cart += z_shift
        # Convert back to fractional with new c
        new_coord = [coord[0], coord[1], z_cart / c_lattice]
        shifted_coords.append(new_coord)
    return shifted_coords

def combine_layers(bot_coords, bot_species, top_coords, top_species):
    """Combine two layers"""
    all_coords = bot_coords + top_coords
    all_species = bot_species + top_species
    
    # Sort by species (Mo, W, S, Se)
    species_order = ['Mo', 'W', 'S', 'Se']
    
    paired = list(zip(all_coords, all_species))
    paired_sorted = sorted(paired, key=lambda x: (species_order.index(x[1]), x))
    
    sorted_coords = [p[0] for p in paired_sorted]
    sorted_species = [p[1] for p in paired_sorted]
    
    # Count species
    species_count = {}
    for sp in sorted_species:
        species_count[sp] = species_count.get(sp, 0) + 1
    
    final_species = [sp for sp in species_order if sp in species_count]
    final_counts = [species_count[sp] for sp in final_species]
    
    return sorted_coords, final_species, final_counts

def write_poscar(filename, structure_name, coords, species, counts, lattice, c_new):
    """Write POSCAR file"""
    with open(filename, 'w') as f:
        f.write(structure_name + "\n")
        f.write("1.0\n")
        
        # Lattice vectors (scale a and b by supercell, use new c)
        f.write("  {0:.10f}  {1:.10f}  {2:.10f}\n".format(
            lattice[0][0], lattice[0][1], lattice[0][2]))
        f.write("  {0:.10f}  {1:.10f}  {2:.10f}\n".format(
            lattice[1][0], lattice[1][1], lattice[1][2]))
        f.write("  {0:.10f}  {1:.10f}  {2:.10f}\n".format(
            0.0, 0.0, c_new))
        
        # Species and counts
        f.write("  " + "  ".join(species) + "\n")
        f.write("  " + "  ".join(map(str, counts)) + "\n")
        
        # Coordinates
        f.write("Direct\n")
        for coord in coords:
            f.write("  {0:.10f}  {1:.10f}  {2:.10f}\n".format(
                coord[0], coord[1], coord[2]))

def write_incar(filename):
    """Write INCAR file"""
    with open(filename, 'w') as f:
        for key, val in incar_template.items():
            if isinstance(val, bool):
                val_str = '.TRUE.' if val else '.FALSE.'
            else:
                val_str = str(val)
            f.write("{0:12s} = {1}\n".format(key, val_str))

def write_kpoints(filename, supercell_size):
    """Write KPOINTS file"""
    # For 5x5 supercell, use 2x2x1 k-mesh
    nx = supercell_size[0]
    k_mesh = max(1, int(9.0 / nx + 0.5))
    
    with open(filename, 'w') as f:
        f.write("Automatic mesh for {0}x{1} supercell\n".format(nx, nx))
        f.write("0\n")
        f.write("Gamma\n")
        f.write("{0} {0} 1\n".format(k_mesh))
        f.write("0 0 0\n")

def build_heterostructure(bot_label, top_label):
    """Build heterostructure from two layer labels"""
    
    print("  Building: {0} // {1}".format(bot_label, top_label))
    
    # Read seed files
    bot_file = structures[bot_label]
    top_file = structures[top_label]
    
    if not os.path.exists(bot_file):
        print("  ERROR: {0} not found".format(bot_file))
        return False
    if not os.path.exists(top_file):
        print("  ERROR: {0} not found".format(top_file))
        return False
    
    bot_struct = read_vasp_file(bot_file)
    top_struct = read_vasp_file(top_file)
    
    print("  Bottom layer: {0} atoms".format(sum(bot_struct['counts'])))
    print("  Top layer: {0} atoms".format(sum(top_struct['counts'])))
    
    # Build supercells
    nx, ny, nz = supercell
    bot_coords, bot_species = make_supercell(bot_struct, nx, ny, nz)
    top_coords, top_species = make_supercell(top_struct, nx, ny, nz)
    
    print("  After supercell: {0} + {1} atoms".format(
        len(bot_coords), len(top_coords)))
    
    # Get original c lattice parameter
    c_orig = bot_struct['lattice'][2][2]
    
    # Position bottom layer at z=0
    bot_coords_shifted = []
    for coord in bot_coords:
        z_cart = coord[2] * c_orig  # Convert to Cartesian
        new_coord = [coord[0], coord[1], z_cart / cell_height_c]
        bot_coords_shifted.append(new_coord)
    
    # Position top layer at z = interlayer_distance
    top_coords_shifted = []
    for coord in top_coords:
        z_cart = coord[2] * c_orig + interlayer_distance
        new_coord = [coord[0], coord[1], z_cart / cell_height_c]
        top_coords_shifted.append(new_coord)
    
    # Combine layers
    all_coords, all_species, all_counts = combine_layers(
        bot_coords_shifted, bot_species,
        top_coords_shifted, top_species
    )
    
    print("  Total atoms: {0}".format(len(all_coords)))
    print("  Species: {0}".format(all_species))
    print("  Counts: {0}".format(all_counts))
    
    # Scale lattice vectors for supercell
    lattice_super = [
        [bot_struct['lattice'][0][0] * nx, 
         bot_struct['lattice'][0][1] * nx, 
         bot_struct['lattice'][0][2]],
        [bot_struct['lattice'][1][0] * ny,
         bot_struct['lattice'][1][1] * ny,
         bot_struct['lattice'][1][2]],
        [0, 0, cell_height_c]
    ]
    
    # Create output directory
    dirname = "{0}__{1}".format(bot_label, top_label)
    if not os.path.exists(dirname):
        os.makedirs(dirname)
    
    # Write files
    structure_name = "{0}__{1} heterostructure".format(bot_label, top_label)
    write_poscar(
        os.path.join(dirname, 'POSCAR'),
        structure_name, all_coords, all_species, all_counts,
        lattice_super, cell_height_c
    )
    write_incar(os.path.join(dirname, 'INCAR'))
    write_kpoints(os.path.join(dirname, 'KPOINTS'), supercell)
    
    print("  SUCCESS: Created {0}/\n".format(dirname))
    return True

# ============================================================
# MAIN
# ============================================================

def main():
    """Main execution"""
    
    print("=" * 60)
    print("Building Janus TMD Heterostructures from Seed Files")
    print("=" * 60)
    print("Supercell: {0}x{1}x{2}".format(*supercell))
    print("Interlayer distance: {0:.2f} A".format(interlayer_distance))
    print("Cell height: {0:.2f} A".format(cell_height_c))
    print()
    
    # Check seed files
    print("Checking seed files...")
    all_exist = True
    for label, filepath in structures.items():
        if os.path.exists(filepath):
            print("  OK: {0}".format(filepath))
        else:
            print("  MISSING: {0}".format(filepath))
            all_exist = False
    print()
    
    if not all_exist:
        print("ERROR: Some seed files are missing!")
        print("Expected files:")
        for label, filepath in structures.items():
            print("  " + filepath)
        sys.exit(1)
    
    # Build heterostructures
    print("Building heterostructures...")
    print()
    
    success_count = 0
    for bot_label, top_label in combinations_to_run:
        if build_heterostructure(bot_label, top_label):
            success_count += 1
    
    print("=" * 60)
    print("DONE! Successfully created {0}/{1} structures".format(
        success_count, len(combinations_to_run)))
    print("=" * 60)
    print()
    print("Next steps:")
    print("1. Copy POTCAR to each directory:")
    print("   for d in MoSSe_*__WSSe_*/; do cp POTCAR \"$d/\"; done")
    print()
    print("2. Verify structures:")
    print("   head -7 MoSSe_Se_up__WSSe_S_up/POSCAR")
    print()
    print("3. Create PBS job scripts and submit")
    print("=" * 60)

if __name__ == '__main__':
    main()
