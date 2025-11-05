#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
fix_broken_poscars_compatible.py
Python 2/3 compatible version

Usage:
  python fix_broken_poscars_compatible.py MoSSe_Se_up__WSSe_S_up
  python fix_broken_poscars_compatible.py MoSSe_S_up__WSSe_Se_up
"""

from __future__ import print_function
import os
import sys

# Z-coordinates configuration (Angstrom)
delta_z = {'S': 1.60, 'Se': 1.65}
d_interlayer = 3.30

# Layer definitions
configs = {
    'MoSSe_Se_up__WSSe_S_up': {
        # Mo layer (z=0): Mo, S_down, Se_up
        # W layer (z=3.3): W, S_down, S_up
        'Mo_z': 0.0,
        'W_z': 3.30,
        'Mo_S_z': -1.60,
        'Mo_Se_z': 1.65,
        'W_S_down_z': 1.70,
        'W_S_up_z': 4.90,
    },
    'MoSSe_S_up__WSSe_Se_up': {
        # Mo layer (z=0): Mo, S_down, S_up
        # W layer (z=3.3): W, S_down, Se_up
        'Mo_z': 0.0,
        'W_z': 3.30,
        'Mo_S_down_z': -1.60,
        'Mo_S_up_z': 1.60,
        'W_S_z': 1.70,
        'W_Se_z': 4.95,
    },
}

def fix_poscar_simple(dirname):
    """Fix POSCAR with correct z-coordinates"""
    
    if dirname not in configs:
        print("Error: Unknown configuration '" + dirname + "'")
        print("Available: " + str(list(configs.keys())))
        return False
    
    config = configs[dirname]
    poscar_path = os.path.join(dirname, 'POSCAR')
    
    if not os.path.exists(poscar_path):
        print("Error: " + poscar_path + " not found")
        return False
    
    # Backup
    backup_path = os.path.join(dirname, 'POSCAR.broken')
    if not os.path.exists(backup_path):
        os.system("cp " + poscar_path + " " + backup_path)
        print("Backed up to " + backup_path)
    
    # Read POSCAR
    with open(poscar_path, 'r') as f:
        lines = f.readlines()
    
    # Parse
    species = lines[5].split()
    counts = list(map(int, lines[6].split()))
    
    print("Species: " + str(species))
    print("Counts: " + str(counts))
    
    # Expected for 5x5 supercell: Mo(25) W(25) S(50) Se(25)
    if species != ['Mo', 'W', 'S', 'Se']:
        print("Warning: Unexpected species order: " + str(species))
        print("Expected: ['Mo', 'W', 'S', 'Se']")
    
    # Get lattice c parameter
    c_val = float(lines[4].split()[2])
    print("Lattice c: {0:.2f} Angstrom".format(c_val))
    
    # Read coordinates
    coord_start = 8
    coords = []
    for i in range(coord_start, coord_start + sum(counts)):
        if i < len(lines):
            parts = lines[i].split()
            if len(parts) >= 3:
                coords.append([float(x) for x in parts[:3]])
    
    print("Read {0} coordinates".format(len(coords)))
    
    # Assign new z-coordinates
    new_coords = []
    atom_idx = 0
    
    for sp, count in zip(species, counts):
        for i in range(count):
            x, y, z_old = coords[atom_idx]
            
            # Determine z based on species and configuration
            if sp == 'Mo':
                z_new = config['Mo_z']
            elif sp == 'W':
                z_new = config['W_z']
            elif sp == 'S':
                # For 5x5 supercell: 25 atoms per chalcogen position
                if dirname == 'MoSSe_Se_up__WSSe_S_up':
                    if i < 25:
                        z_new = config['Mo_S_z']
                    elif i < 37:
                        z_new = config['W_S_down_z']
                    else:
                        z_new = config['W_S_up_z']
                else:  # MoSSe_S_up__WSSe_Se_up
                    if i < 12:
                        z_new = config['Mo_S_down_z']
                    elif i < 25:
                        z_new = config['Mo_S_up_z']
                    else:
                        z_new = config['W_S_z']
            elif sp == 'Se':
                if dirname == 'MoSSe_Se_up__WSSe_S_up':
                    z_new = config['Mo_Se_z']
                else:  # MoSSe_S_up__WSSe_Se_up
                    z_new = config['W_Se_z']
            else:
                z_new = z_old * c_val
            
            # Convert to fractional
            z_frac = z_new / c_val
            new_coords.append([x, y, z_frac])
            atom_idx += 1
    
    # Write fixed POSCAR
    with open(poscar_path, 'w') as f:
        # Write header (unchanged)
        for i in range(coord_start):
            f.write(lines[i])
        
        # Write new coordinates
        for coord in new_coords:
            f.write("  {0:.10f}  {1:.10f}  {2:.10f}\n".format(
                coord[0], coord[1], coord[2]))
    
    print("Fixed POSCAR written to " + poscar_path)
    
    # Verify
    z_values = sorted(set("{0:.6f}".format(c[2]) for c in new_coords))
    print("Unique z-coordinates: {0}".format(len(z_values)))
    if len(z_values) <= 10:
        print("Z-values: " + str(z_values))
    else:
        print("Z-values: " + str(z_values[:10]) + "...")
    
    if len(z_values) < 4:
        print("WARNING: Less than 4 unique z-values. Check structure!")
        return False
    
    return True

def main():
    """Main function"""
    
    if len(sys.argv) < 2:
        print("Usage: python fix_broken_poscars_compatible.py <directory>")
        print("")
        print("Available directories:")
        for name in configs.keys():
            print("  " + name)
        sys.exit(1)
    
    dirname = sys.argv[1]
    
    print("=" * 60)
    print("Fixing POSCAR in: " + dirname)
    print("=" * 60)
    print()
    
    success = fix_poscar_simple(dirname)
    
    print()
    print("=" * 60)
    if success:
        print("SUCCESS!")
        print()
        print("Next steps:")
        print("  1. Check structure: head -30 " + dirname + "/POSCAR")
        print("  2. Verify z-coords: awk 'NR>8 {print $3}' " + dirname + "/POSCAR | sort -u")
        print("  3. If OK, resubmit: cd " + dirname + " && qsub vasp.job")
    else:
        print("FAILED - Please check errors above")
    print("=" * 60)

if __name__ == '__main__':
    main()
