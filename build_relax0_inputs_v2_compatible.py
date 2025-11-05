#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
build_relax0_inputs_v2_fixed_compatible.py
Python 2/3 compatible version - no f-strings
"""

from __future__ import print_function
import os
import sys

# ========== Configuration ==========
delta_z = {
    'S': 1.60,   # Angstrom
    'Se': 1.65,  # Angstrom
}

d_interlayer = 3.30  # Angstrom, Mo-W distance

# More robust INCAR template
incar_template = {
    'PREC': 'Normal',
    'ENCUT': 400,
    'ALGO': 'Normal',
    'NELM': 150,
    'EDIFF': 1e-5,
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
    'NWRITE': 2,
    'NCORE': 4,
}

# All stacking configurations
stacking_configs = [
    ('MoS2_up', 'WS2_up'),
    ('MoS2_up', 'WS2_down'),
    ('MoSSe_S_up', 'WS2_up'),
    ('MoSSe_S_up', 'WS2_down'),
    ('MoSSe_Se_up', 'WS2_up'),
    ('MoSSe_Se_up', 'WS2_down'),
    ('MoS2_up', 'WSSe_S_up'),
    ('MoS2_up', 'WSSe_Se_up'),
    ('MoSSe_S_up', 'WSSe_S_up'),
    ('MoSSe_S_up', 'WSSe_Se_up'),
    ('MoSSe_Se_up', 'WSSe_S_up'),
    ('MoSSe_Se_up', 'WSSe_Se_up'),
]

def parse_layer_label(label):
    """Parse layer configuration label"""
    if label.startswith('MoS2'):
        metal, X_base = 'Mo', 'S'
        X_top = 'S'
    elif label.startswith('MoSSe'):
        metal, X_base = 'Mo', 'S'
        parts = label.split('_')
        X_top = parts[1]  # 'S' or 'Se'
    elif label.startswith('WS2'):
        metal, X_base = 'W', 'S'
        X_top = 'S'
    elif label.startswith('WSSe'):
        metal, X_base = 'W', 'S'
        parts = label.split('_')
        X_top = parts[1]
    else:
        raise ValueError("Unknown label: " + label)
    
    direction = label.split('_')[-1]  # 'up' or 'down'
    return metal, X_base, X_top, direction

def build_layer_coords(metal, X_base, X_top, direction, z_metal, a):
    """Build atomic coordinates for one layer"""
    coords = []
    types = []
    
    # Metal atom
    coords.append([a/3.0, a/3.0, z_metal])
    types.append(metal)
    
    # Calculate chalcogen positions
    if direction == 'up':
        z_down = z_metal - delta_z[X_base]
        z_up = z_metal + delta_z[X_top]
        X_down, X_up = X_base, X_top
    else:  # down
        z_down = z_metal - delta_z[X_top]
        z_up = z_metal + delta_z[X_base]
        X_down, X_up = X_top, X_base
    
    # Bottom chalcogen
    coords.append([0.0, 0.0, z_down])
    types.append(X_down)
    
    # Top chalcogen
    coords.append([2*a/3.0, 2*a/3.0, z_up])
    types.append(X_up)
    
    return coords, types

def check_structure_validity(coords, types, min_distance=1.5):
    """Check if structure has reasonable atomic distances"""
    n_atoms = len(coords)
    
    for i in range(n_atoms):
        for j in range(i+1, n_atoms):
            dz = abs(coords[i][2] - coords[j][2])
            
            if dz < min_distance:
                msg = "WARNING: Atoms {0}({1}) and {2}({3}) too close: dz={4:.3f} A".format(
                    i, types[i], j, types[j], dz)
                return False, msg
    
    return True, "Structure OK"

def write_poscar(filepath, layer1_label, layer2_label, a, c_lattice):
    """Generate POSCAR with structure validation"""
    # Parse layers
    M1, X1_base, X1_top, dir1 = parse_layer_label(layer1_label)
    M2, X2_base, X2_top, dir2 = parse_layer_label(layer2_label)
    
    # Build coordinates
    z_Mo = 0.0
    z_W = d_interlayer
    coords1, types1 = build_layer_coords(M1, X1_base, X1_top, dir1, z_Mo, a)
    coords2, types2 = build_layer_coords(M2, X2_base, X2_top, dir2, z_W, a)
    
    all_coords = coords1 + coords2
    all_types = types1 + types2
    
    # Validate structure
    is_valid, msg = check_structure_validity(all_coords, all_types)
    if not is_valid:
        print("  " + msg)
    
    # Sort by species
    species_order = ['Mo', 'W', 'S', 'Se']
    sorted_indices = sorted(range(len(all_types)), 
                          key=lambda i: (species_order.index(all_types[i]), i))
    
    sorted_coords = [all_coords[i] for i in sorted_indices]
    sorted_types = [all_types[i] for i in sorted_indices]
    
    # Count atoms
    species_count = {}
    for t in sorted_types:
        species_count[t] = species_count.get(t, 0) + 1
    
    species_list = [sp for sp in species_order if sp in species_count]
    counts = [species_count[sp] for sp in species_list]
    
    # Write POSCAR
    with open(filepath, 'w') as f:
        f.write(layer1_label + "__" + layer2_label + "\n")
        f.write("1.0\n")
        f.write("  {0:.10f}  {1:.10f}  {2:.10f}\n".format(a, 0.0, 0.0))
        f.write("  {0:.10f}  {1:.10f}  {2:.10f}\n".format(-a/2, a*0.866025404/2, 0.0))
        f.write("  {0:.10f}  {1:.10f}  {2:.10f}\n".format(0.0, 0.0, c_lattice))
        f.write("  " + "  ".join(species_list) + "\n")
        f.write("  " + "  ".join(map(str, counts)) + "\n")
        f.write("Direct\n")
        
        for coord in sorted_coords:
            # Convert to fractional coordinates
            frac_x = coord[0] / a
            frac_y = coord[1] / (a * 0.866025404)
            frac_z = coord[2] / c_lattice
            f.write("  {0:.10f}  {1:.10f}  {2:.10f}\n".format(frac_x, frac_y, frac_z))

def write_incar(filepath):
    """Write INCAR file"""
    with open(filepath, 'w') as f:
        for key, val in incar_template.items():
            if isinstance(val, bool):
                val_str = '.TRUE.' if val else '.FALSE.'
            else:
                val_str = str(val)
            f.write("{0:12s} = {1}\n".format(key, val_str))

def write_kpoints(filepath):
    """Write KPOINTS file"""
    with open(filepath, 'w') as f:
        f.write("Automatic mesh\n")
        f.write("0\n")
        f.write("Gamma\n")
        f.write("9 9 1\n")
        f.write("0 0 0\n")

def main():
    """Main execution"""
    # Lattice parameters
    a = 3.18  # in-plane lattice constant (Angstrom)
    c_lattice = 30.0  # large c to avoid interaction
    
    print("=" * 60)
    print("Building VASP input files (IMPROVED VERSION)")
    print("=" * 60)
    print("Chalcogen z-spacings: S={0:.2f} A, Se={1:.2f} A".format(
        delta_z['S'], delta_z['Se']))
    print("Interlayer distance: {0:.2f} A".format(d_interlayer))
    print()
    
    # Create directories and files
    for layer1, layer2 in stacking_configs:
        dirname = layer1 + "__" + layer2
        
        if not os.path.exists(dirname):
            os.makedirs(dirname)
        
        print("Creating: " + dirname)
        
        # Write files
        write_poscar(os.path.join(dirname, 'POSCAR'), layer1, layer2, a, c_lattice)
        write_incar(os.path.join(dirname, 'INCAR'))
        write_kpoints(os.path.join(dirname, 'KPOINTS'))
    
    print()
    print("=" * 60)
    print("DONE! Remember to:")
    print("1. Copy POTCAR to each directory")
    print("2. Check structures: grep 'WARNING' */POSCAR")
    print("3. Consider 2-stage relaxation for problematic cases")
    print("=" * 60)

if __name__ == '__main__':
    main()
