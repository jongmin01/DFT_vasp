#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
generate_potcar_from_Z.py
Generate POTCAR from compressed POTCAR.Z files
Python 2/3 compatible
"""

from __future__ import print_function
import os
import sys
import subprocess

# ============================================================
# CONFIGURATION
# ============================================================

# POTCAR base directory
POTCAR_BASE = "/opt/soft/vasp-pot/ase/potpaw_PBE"

# Elements in order (must match POSCAR)
elements = ['Mo_pv', 'W_pv', 'S', 'Se']

# ============================================================
# FUNCTIONS
# ============================================================

def run_command(cmd):
    """Run shell command and return output"""
    try:
        output = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
        if isinstance(output, bytes):  # Python 3
            output = output.decode('utf-8')
        return output.strip()
    except subprocess.CalledProcessError as e:
        return None

def check_potcar_exists(element, base_dir):
    """Check if POTCAR.Z file exists for element"""
    potcar_path = os.path.join(base_dir, element, 'POTCAR.Z')
    return os.path.exists(potcar_path)

def generate_potcar(output_file='POTCAR'):
    """Generate combined POTCAR file from .Z files"""
    
    print("=" * 60)
    print("POTCAR Generator (from .Z files)")
    print("=" * 60)
    print("Base directory: " + POTCAR_BASE)
    print("Elements: " + str(elements))
    print()
    
    # Check if base directory exists
    if not os.path.exists(POTCAR_BASE):
        print("ERROR: POTCAR base directory not found:")
        print("  " + POTCAR_BASE)
        print()
        print("Please check the path or modify POTCAR_BASE in this script")
        sys.exit(1)
    
    # Check all POTCAR.Z files
    print("Checking POTCAR.Z files...")
    all_exist = True
    potcar_files = []
    
    for elem in elements:
        potcar_path = os.path.join(POTCAR_BASE, elem, 'POTCAR.Z')
        if os.path.exists(potcar_path):
            print("  OK: {0:10s} -> {1}".format(elem, potcar_path))
            potcar_files.append(potcar_path)
        else:
            print("  MISSING: " + potcar_path)
            all_exist = False
    
    print()
    
    if not all_exist:
        print("ERROR: Some POTCAR.Z files are missing!")
        sys.exit(1)
    
    # Check if zcat is available
    zcat_check = run_command("which zcat")
    if not zcat_check:
        print("ERROR: 'zcat' command not found!")
        print("This script requires zcat to decompress .Z files")
        sys.exit(1)
    
    print("Generating POTCAR...")
    
    # Remove existing POTCAR if present
    if os.path.exists(output_file):
        os.remove(output_file)
        print("  Removed existing " + output_file)
    
    # Combine POTCAR files using zcat
    for i, potcar_path in enumerate(potcar_files):
        print("  Adding {0}/{1}: {2}...".format(
            i+1, len(potcar_files), elements[i]))
        
        cmd = "zcat {0} >> {1}".format(potcar_path, output_file)
        result = run_command(cmd)
        
        if result is None:
            # Command succeeded (no output expected)
            pass
        else:
            print("    Warning: " + str(result))
    
    print()
    print("POTCAR generated: " + output_file)
    print()
    
    # Verify
    print("Verification (TITEL lines):")
    if os.path.exists(output_file):
        with open(output_file, 'r') as f:
            for line in f:
                if 'TITEL' in line:
                    print("  " + line.strip())
        
        # File size
        size = os.path.getsize(output_file)
        print()
        print("POTCAR size: {0} bytes (~{1} KB)".format(size, size // 1024))
    else:
        print("  ERROR: POTCAR was not created!")
        sys.exit(1)
    
    print()
    print("=" * 60)
    print("SUCCESS!")
    print("=" * 60)
    print()
    print("Next steps:")
    print("1. Copy to all calculation directories:")
    print("   for d in MoSSe_*__WSSe_*/; do cp POTCAR \"$d/\"; done")
    print()
    print("2. Verify in each directory:")
    print("   ls -lh MoSSe_*/POTCAR")

def main():
    """Main function"""
    
    if len(sys.argv) > 1:
        if sys.argv[1] in ['-h', '--help']:
            print("Usage: python generate_potcar_from_Z.py [output_file]")
            print()
            print("Default output: POTCAR")
            print()
            print("Edit POTCAR_BASE in script to set your POTCAR location")
            sys.exit(0)
        output_file = sys.argv[1]
    else:
        output_file = 'POTCAR'
    
    generate_potcar(output_file)

if __name__ == '__main__':
    main()
