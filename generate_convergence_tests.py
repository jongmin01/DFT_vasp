#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
generate_convergence_tests.py
Generate convergence test calculations for k-points and ENCUT
"""

from __future__ import print_function
import os
import sys

# ============================================================
# CONVERGENCE TEST SETTINGS
# ============================================================

# K-point meshes to test
kpoint_meshes = [
    (2, 2, 1),
    (3, 3, 1),
    (4, 4, 1),
    (5, 5, 1),
]

# ENCUT values to test (eV)
encut_values = [400, 450, 500, 520, 550]

# Base directory for tests
test_base_dir = "convergence_tests"

# Reference structure to use
reference_structure = "MoSSe_Se_up__WSSe_S_up"

# ============================================================
# INCAR TEMPLATE FOR CONVERGENCE TESTS
# ============================================================

# Single-point calculation (no relaxation)
incar_convergence_base = {
    'SYSTEM': 'Convergence test',
    'PREC': 'Accurate',
    'ENCUT': 520,  # Will be overridden
    'EDIFF': 1e-7,  # Tight for accurate energies
    'ALGO': 'Normal',
    'NELM': 200,
    
    # No relaxation
    'IBRION': -1,
    'NSW': 0,
    
    # Electronic settings
    'ISPIN': 2,
    'MAGMOM': '25*2.0 25*2.0 100*0.0',
    'IVDW': 12,
    'ISMEAR': 0,
    'SIGMA': 0.05,
    
    # Accuracy
    'LASPH': '.TRUE.',
    'LMAXMIX': 4,
    'ADDGRID': '.TRUE.',
    'LREAL': '.FALSE.',
    
    # Output
    'LORBIT': 11,
    'LWAVE': '.FALSE.',
    'LCHARG': '.FALSE.',
    'NCORE': 4,
}

# ============================================================
# FUNCTIONS
# ============================================================

def write_incar(filename, incar_dict):
    """Write INCAR file"""
    with open(filename, 'w') as f:
        for key, val in incar_dict.items():
            if isinstance(val, bool):
                val_str = '.TRUE.' if val else '.FALSE.'
            elif isinstance(val, float):
                if abs(val) < 0.01 and val != 0:
                    val_str = "{0:.2E}".format(val)
                else:
                    val_str = str(val)
            else:
                val_str = str(val)
            f.write("{0:12s} = {1}\n".format(key, val_str))

def write_kpoints(filename, comment, mesh):
    """Write KPOINTS file"""
    with open(filename, 'w') as f:
        f.write(comment + "\n")
        f.write("0\n")
        f.write("Gamma\n")
        f.write("{0} {1} {2}\n".format(mesh[0], mesh[1], mesh[2]))
        f.write("0 0 0\n")

def generate_kpoint_tests():
    """Generate k-point convergence tests"""
    
    test_dir = os.path.join(test_base_dir, "kpoint_convergence")
    
    if not os.path.exists(test_dir):
        os.makedirs(test_dir)
    
    print("Generating k-point convergence tests...")
    print("  Location: {0}".format(test_dir))
    print()
    
    for kx, ky, kz in kpoint_meshes:
        # Create directory
        dirname = "ktest_{0}x{1}x{2}".format(kx, ky, kz)
        full_path = os.path.join(test_dir, dirname)
        
        if not os.path.exists(full_path):
            os.makedirs(full_path)
        
        print("  Creating: {0}".format(dirname))
        
        # Copy POSCAR and POTCAR from reference
        for fname in ['POSCAR', 'POTCAR']:
            src = os.path.join(reference_structure, fname)
            dst = os.path.join(full_path, fname)
            if os.path.exists(src):
                os.system("cp {0} {1}".format(src, dst))
        
        # Write INCAR
        incar = incar_convergence_base.copy()
        incar['SYSTEM'] = "K-point test {0}x{1}x{2}".format(kx, ky, kz)
        write_incar(os.path.join(full_path, 'INCAR'), incar)
        
        # Write KPOINTS
        comment = "K-point convergence test {0}x{1}x{2}".format(kx, ky, kz)
        write_kpoints(os.path.join(full_path, 'KPOINTS'), comment, (kx, ky, kz))
    
    print()
    print("  Created {0} k-point tests".format(len(kpoint_meshes)))
    print()

def generate_encut_tests():
    """Generate ENCUT convergence tests"""
    
    test_dir = os.path.join(test_base_dir, "encut_convergence")
    
    if not os.path.exists(test_dir):
        os.makedirs(test_dir)
    
    print("Generating ENCUT convergence tests...")
    print("  Location: {0}".format(test_dir))
    print()
    
    for encut in encut_values:
        # Create directory
        dirname = "encut_{0}".format(encut)
        full_path = os.path.join(test_dir, dirname)
        
        if not os.path.exists(full_path):
            os.makedirs(full_path)
        
        print("  Creating: {0}".format(dirname))
        
        # Copy POSCAR and POTCAR from reference
        for fname in ['POSCAR', 'POTCAR']:
            src = os.path.join(reference_structure, fname)
            dst = os.path.join(full_path, fname)
            if os.path.exists(src):
                os.system("cp {0} {1}".format(src, dst))
        
        # Write INCAR
        incar = incar_convergence_base.copy()
        incar['SYSTEM'] = "ENCUT test {0} eV".format(encut)
        incar['ENCUT'] = encut
        write_incar(os.path.join(full_path, 'INCAR'), incar)
        
        # Write KPOINTS (use 3x3x1 for all ENCUT tests)
        comment = "ENCUT convergence test {0} eV".format(encut)
        write_kpoints(os.path.join(full_path, 'KPOINTS'), comment, (3, 3, 1))
    
    print()
    print("  Created {0} ENCUT tests".format(len(encut_values)))
    print()

def generate_analysis_script():
    """Generate script to analyze convergence results"""
    
    script_content = '''#!/bin/bash
# analyze_convergence.sh
# Extract and compare energies from convergence tests

echo "======================================================"
echo "Convergence Test Analysis"
echo "======================================================"
echo ""

# K-point convergence
if [ -d "convergence_tests/kpoint_convergence" ]; then
    echo "K-point Convergence:"
    echo "  Mesh      Energy (eV)    E-E_ref (meV/atom)"
    echo "  ----      -----------    ------------------"
    
    cd convergence_tests/kpoint_convergence
    
    # Get reference energy (finest mesh)
    ref_dir=$(ls -d ktest_* | tail -1)
    ref_energy=$(grep "free  energy" $ref_dir/OUTCAR 2>/dev/null | tail -1 | awk '{print $5}')
    
    for d in ktest_*/; do
        mesh=$(echo $d | sed 's/ktest_//;s/\///')
        
        if [ -f "$d/OUTCAR" ]; then
            energy=$(grep "free  energy" $d/OUTCAR 2>/dev/null | tail -1 | awk '{print $5}')
            
            if [ -n "$energy" ] && [ -n "$ref_energy" ]; then
                # Calculate difference per atom (150 atoms)
                diff=$(echo "scale=6; ($energy - $ref_energy) * 1000 / 150" | bc)
                printf "  %-10s %-15s %10.3f\\n" "$mesh" "$energy" "$diff"
            else
                printf "  %-10s %-15s %10s\\n" "$mesh" "Not finished" "---"
            fi
        else
            printf "  %-10s %-15s %10s\\n" "$mesh" "Not run" "---"
        fi
    done
    
    cd ../..
    echo ""
fi

# ENCUT convergence
if [ -d "convergence_tests/encut_convergence" ]; then
    echo "ENCUT Convergence:"
    echo "  ENCUT(eV)  Energy (eV)    E-E_ref (meV/atom)"
    echo "  ---------  -----------    ------------------"
    
    cd convergence_tests/encut_convergence
    
    # Get reference energy (highest ENCUT)
    ref_dir=$(ls -d encut_* | sort -t_ -k2 -n | tail -1)
    ref_energy=$(grep "free  energy" $ref_dir/OUTCAR 2>/dev/null | tail -1 | awk '{print $5}')
    
    for d in $(ls -d encut_* | sort -t_ -k2 -n); do
        encut=$(echo $d | sed 's/encut_//;s/\///')
        
        if [ -f "$d/OUTCAR" ]; then
            energy=$(grep "free  energy" $d/OUTCAR 2>/dev/null | tail -1 | awk '{print $5}')
            
            if [ -n "$energy" ] && [ -n "$ref_energy" ]; then
                # Calculate difference per atom
                diff=$(echo "scale=6; ($energy - $ref_energy) * 1000 / 150" | bc)
                printf "  %-10s %-15s %10.3f\\n" "$encut" "$energy" "$diff"
            else
                printf "  %-10s %-15s %10s\\n" "$encut" "Not finished" "---"
            fi
        else
            printf "  %-10s %-15s %10s\\n" "$encut" "Not run" "---"
        fi
    done
    
    cd ../..
    echo ""
fi

echo "======================================================"
echo "Convergence Criteria:"
echo "  K-points: Difference < 1 meV/atom"
echo "  ENCUT:    Difference < 1 meV/atom"
echo "======================================================"
'''
    
    with open('analyze_convergence.sh', 'w') as f:
        f.write(script_content)
    
    os.chmod('analyze_convergence.sh', 0o755)
    print("Created: analyze_convergence.sh")

def generate_submit_script():
    """Generate script to submit convergence tests"""
    
    script_content = '''#!/bin/bash
# submit_convergence_tests.sh
# Submit all convergence test jobs

echo "Submitting convergence tests..."
echo ""

# K-point tests (quick, use 1 node)
if [ -d "convergence_tests/kpoint_convergence" ]; then
    echo "K-point tests:"
    cd convergence_tests/kpoint_convergence
    for d in ktest_*/; do
        if [ -f "$d/run.pbs" ]; then
            echo "  Submitting: $d"
            cd "$d"
            qsub run.pbs
            cd ..
        fi
    done
    cd ../..
    echo ""
fi

# ENCUT tests (quick, use 1 node)
if [ -d "convergence_tests/encut_convergence" ]; then
    echo "ENCUT tests:"
    cd convergence_tests/encut_convergence
    for d in encut_*/; do
        if [ -f "$d/run.pbs" ]; then
            echo "  Submitting: $d"
            cd "$d"
            qsub run.pbs
            cd ..
        fi
    done
    cd ../..
    echo ""
fi

echo "All convergence tests submitted"
echo "Check status: qstat -u \\$USER"
echo "Analyze results: ./analyze_convergence.sh"
'''
    
    with open('submit_convergence_tests.sh', 'w') as f:
        f.write(script_content)
    
    os.chmod('submit_convergence_tests.sh', 0o755)
    print("Created: submit_convergence_tests.sh")

# ============================================================
# MAIN
# ============================================================

def main():
    """Main function"""
    
    print("=" * 60)
    print("Convergence Test Generator")
    print("=" * 60)
    print()
    
    # Check reference structure exists
    if not os.path.isdir(reference_structure):
        print("ERROR: Reference structure not found: " + reference_structure)
        print()
        print("Available directories:")
        for d in os.listdir('.'):
            if os.path.isdir(d):
                print("  " + d)
        sys.exit(1)
    
    print("Reference structure: {0}".format(reference_structure))
    print()
    
    # Generate tests
    generate_kpoint_tests()
    generate_encut_tests()
    
    # Generate analysis scripts
    generate_analysis_script()
    generate_submit_script()
    
    print()
    print("=" * 60)
    print("Convergence tests generated!")
    print("=" * 60)
    print()
    print("Structure:")
    print("  convergence_tests/")
    print("    kpoint_convergence/")
    print("      ktest_2x2x1/")
    print("      ktest_3x3x1/")
    print("      ktest_4x4x1/")
    print("      ktest_5x5x1/")
    print("    encut_convergence/")
    print("      encut_400/")
    print("      encut_450/")
    print("      encut_500/")
    print("      encut_520/")
    print("      encut_550/")
    print()
    print("Next steps:")
    print("  1. Generate PBS scripts for tests:")
    print("     cd convergence_tests/kpoint_convergence")
    print("     for d in */; do")
    print("       # Create run.pbs in each directory")
    print("     done")
    print()
    print("  2. Or manually copy from main calculations:")
    print("     for d in convergence_tests/*/*/; do")
    print("       cp {0}/run.pbs \"$d/\"".format(reference_structure))
    print("       # Edit walltime to ~2 hours (quick tests)")
    print("     done")
    print()
    print("  3. Submit tests:")
    print("     ./submit_convergence_tests.sh")
    print()
    print("  4. Analyze results:")
    print("     ./analyze_convergence.sh")
    print()
    print("=" * 60)

if __name__ == '__main__':
    main()
