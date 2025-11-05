#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
build_relax0_inputs_v2.py
----------------------------------------------------
Generate 0-degree Janus heterostructure (MoSSe–WSSe)
and prepare VASP relaxation input files.

Improvements:
- Proper MAGMOM settings
- Gamma-centered KPOINTS
- Optimized k-mesh for supercell
- ISIF=3 for full relaxation
- Proper NPAR calculation
- All combinations of 4 structures

Author: Updated version based on discussions
Date: 2025-11-01
----------------------------------------------------
"""

import os
import numpy as np
from ase.io import read, write
from ase.build import stack
from collections import Counter

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
interlayer_distance = 6.8           # Interlayer distance in Å (initial guess)
cell_height_c = 30.0                # Total cell height in Å

# POTCAR information
potcar_path = "./POTCAR"            # POTCAR file path
# Expected POTCAR order: Mo_sv W_sv S Se (verify with: grep TITEL POTCAR)

# Output base directory
output_base = "./relax_0deg_structures"

# PBS and VASP parameters
pbs_account = "cnm84150"
pbs_queue = "batch"
pbs_nodes = 3
pbs_ppn = 16
pbs_gen = "gen6"               # gen6 has 192GB RAM (needed for 150 atoms)
walltime = "72:00:00"
vasp_module = "vasp5"          # or "vasp/5.4.4-intel"
email = "jbaek27@uic.edu"
notify_mode = "abe"

# Calculation settings
calculation_type = "relax"     # "relax" or "test"
k_mesh_density = "normal"      # "coarse", "normal", "fine"

# ============================================================
# COMBINATIONS TO CALCULATE
# ============================================================

# Define which combinations to calculate
# Format: (bottom_layer, top_layer)
combinations_to_run = [
    ('MoSSe_S_up', 'WSSe_S_up'),      # Example 1
    ('MoSSe_S_up', 'WSSe_Se_up'),    # Example 2
    # Add more combinations as needed:
    ('MoSSe_Se_up', 'WSSe_S_up'),
    ('MoSSe_Se_up', 'WSSe_Se_up'),
    # ... (total 16 possible combinations)
]

# ============================================================
# HELPER FUNCTIONS
# ============================================================

def calculate_optimal_kpoints(supercell_size, density='normal'):
    """
    Calculate optimal k-point mesh for supercell
    
    Parameters:
    -----------
    supercell_size : int
        Supercell expansion factor (e.g., 5 for 5×5)
    density : str
        'coarse', 'normal', or 'fine'
    
    Returns:
    --------
    k_mesh : int
        k-point mesh size (for n×n×1)
    """
    
    # Base k-mesh for primitive cell
    k_base_dict = {
        'coarse': 12,
        'normal': 15,
        'fine': 18
    }
    
    k_base = k_base_dict.get(density, 15)
    k_mesh = max(2, k_base // supercell_size)
    
    return k_mesh


def calculate_npar(n_cores):
    """
    Calculate optimal NPAR for given number of cores
    Rule: NPAR should divide n_cores, aim for ~8 bands per group
    """
    
    # Try values: 4, 6, 8, 12, 16
    candidates = [4, 6, 8, 12, 16]
    
    for npar in candidates:
        if n_cores % npar == 0:
            return npar
    
    # Fallback: largest divisor
    for npar in range(8, 2, -1):
        if n_cores % npar == 0:
            return npar
    
    return 4  # Safe default


def analyze_structure(atoms):
    """
    Analyze atomic structure and return element counts and order
    """
    
    symbols = atoms.get_chemical_symbols()
    element_counts = Counter(symbols)
    
    # Get unique elements in order of appearance
    element_order = []
    seen = set()
    for sym in symbols:
        if sym not in seen:
            element_order.append(sym)
            seen.add(sym)
    
    return element_counts, element_order


def generate_magmom(element_counts, element_order):
    """
    Generate MAGMOM string based on element types
    
    TMDs: Mo and W have small magnetic moments
          S and Se are non-magnetic
    """
    
    magmom_values = {
        'Mo': 0.5,
        'W': 0.5,
        'S': 0.0,
        'Se': 0.0
    }
    
    magmom_parts = []
    for element in element_order:
        count = element_counts[element]
        mag = magmom_values.get(element, 0.0)
        magmom_parts.append(f"{count}*{mag}")
    
    return " ".join(magmom_parts)


# ============================================================
# INPUT FILE GENERATORS
# ============================================================

def generate_incar(atoms, n_cores, calculation_type='relax'):
    """
    Generate INCAR file with proper settings
    
    Parameters:
    -----------
    atoms : ASE Atoms object
        Structure to get element information
    n_cores : int
        Total number of cores
    calculation_type : str
        'relax' or 'test'
    """
    
    # Analyze structure
    element_counts, element_order = analyze_structure(atoms)
    magmom_string = generate_magmom(element_counts, element_order)
    
    # Calculate NPAR
    npar = calculate_npar(n_cores)
    
    # Adjust parameters based on calculation type
    if calculation_type == 'test':
        nsw = 50
        ediffg = -0.02
    else:  # relax
        nsw = 300
        ediffg = -0.015
    
    incar = f"""SYSTEM  = Janus TMD Heterostructure 0deg relax

# === Electronic Structure ===
ENCUT   = 550          # Cutoff energy (eV)
PREC    = Accurate     # Precision level
EDIFF   = 1E-6         # Electronic convergence
ALGO    = Normal       # Electronic minimization algorithm
NELM    = 150          # Max electronic steps

# === Ionic Relaxation ===
IBRION  = 2            # CG algorithm
NSW     = {nsw}            # Max ionic steps
EDIFFG  = {ediffg}        # Force convergence (eV/Angstrom)
POTIM   = 0.5          # Step size scaling

# === Cell Relaxation ===
ISIF    = 3            # Relax ions + cell shape + volume

# === Spin Polarization ===
ISPIN   = 2            # Spin-polarized calculation
MAGMOM  = {magmom_string}

# === vdW Correction (CRITICAL for layered materials!) ===
IVDW    = 12           # DFT-D3(BJ) method

# === Accuracy Improvements ===
LASPH   = .TRUE.       # Non-spherical contributions (important for d-electrons)
LMAXMIX = 4            # For d-electrons
ADDGRID = .TRUE.       # Finer grid for augmentation charges

# === k-space Integration ===
ISMEAR  = 0            # Gaussian smearing (good for semiconductors)
SIGMA   = 0.05         # Smearing width (eV)

# === Parallelization ===
NPAR    = {npar}           # Parallelization over bands
LREAL   = Auto         # Real-space projection

# === Output Control ===
LORBIT  = 11           # Write DOSCAR and projected DOS
LWAVE   = .FALSE.      # Don't write WAVECAR (save space)
LCHARG  = .FALSE.      # Don't write CHGCAR (save space)

# === Additional Settings ===
NEDOS   = 2000         # DOS points
"""
    
    return incar


def generate_kpoints(supercell_size, density='normal', calc_type='relax'):
    """
    Generate KPOINTS file
    Always use Gamma-centered for 2D systems!
    
    Parameters:
    -----------
    supercell_size : int
        Supercell expansion factor
    density : str
        'coarse', 'normal', or 'fine'
    calc_type : str
        'relax', 'scf', or 'dos'
    """
    
    if calc_type == 'relax':
        k_mesh = calculate_optimal_kpoints(supercell_size, density)
    elif calc_type == 'scf':
        k_mesh = calculate_optimal_kpoints(supercell_size, 'fine')
        k_mesh = max(k_mesh, 5)
    elif calc_type == 'dos':
        k_mesh = calculate_optimal_kpoints(supercell_size, 'fine')
        k_mesh = max(k_mesh + 2, 7)
    else:
        k_mesh = 3
    
    kpoints = f"""Automatic mesh - {calc_type.upper()}
0
Gamma
{k_mesh} {k_mesh} 1
0 0 0
"""
    
    return kpoints, k_mesh


def generate_pbs_script(
    job_name,
    account,
    queue,
    nodes,
    ppn,
    gen,
    walltime,
    vasp_module,
    email,
    notify_mode,
    work_dir
):
    """
    Generate PBS job script for Carbon cluster
    """
    
    n_cores = nodes * ppn
    
    pbs_script = f"""#!/bin/bash
#PBS -N {job_name}
#PBS -A {account}
#PBS -q {queue}
#PBS -l nodes={nodes}:ppn={ppn}:{gen}
#PBS -l walltime={walltime}
#PBS -l naccesspolicy=SINGLEJOB -n
#PBS -j oe
#PBS -m {notify_mode}
#PBS -M {email}

# ============================================================
# PBS Job Script for VASP Calculation
# Generated automatically
# ============================================================

echo "========================================"
echo "Job started on $(hostname) at $(date)"
echo "Job ID: $PBS_JOBID"
echo "Job Name: $PBS_JOBNAME"
echo "Working directory: $PBS_O_WORKDIR"
echo "========================================"

# Load modules
module load {vasp_module}
module list

# Change to working directory
cd $PBS_O_WORKDIR
echo "Running in: $(pwd)"

# Print node information
echo "Nodes allocated:"
cat $PBS_NODEFILE
echo ""
echo "Number of cores: {n_cores}"
echo ""

# Print input file info
echo "========================================"
echo "INCAR settings:"
grep -E "ISIF|ISPIN|IVDW|ENCUT|NSW" INCAR
echo ""
echo "KPOINTS:"
head -4 KPOINTS
echo "========================================"

# Run VASP
echo "Starting VASP calculation at $(date)"
time mpirun -np {n_cores} vasp_std > vasp.log 2>&1

# Check completion
if [ -f CONTCAR ]; then
    echo "Calculation completed at $(date)"
    
    # Quick check
    if grep -q "reached required accuracy" OUTCAR; then
        echo "✓ Convergence achieved!"
    else
        echo "⚠ Check convergence in OUTCAR"
    fi
    
    # Extract final energy
    if [ -f OSZICAR ]; then
        echo ""
        echo "Final energies:"
        tail -5 OSZICAR
    fi
else
    echo "⚠ Warning: CONTCAR not found. Calculation may have failed."
fi

echo "========================================"
echo "Job finished at $(date)"
echo "========================================"
"""
    
    return pbs_script


# ============================================================
# MAIN WORKFLOW
# ============================================================

def build_heterostructure(bot_file, top_file, supercell, distance, cell_c):
    """
    Build heterostructure from two layer files
    """
    
    print(f"  Reading bottom: {os.path.basename(bot_file)}")
    print(f"  Reading top: {os.path.basename(top_file)}")
    
    bot = read(bot_file)
    top = read(top_file)
    
    # Build supercell
    bot_super = bot.repeat(supercell)
    top_super = top.repeat(supercell)
    
    print(f"  Supercell: {supercell[0]}×{supercell[1]}×{supercell[2]}")
    print(f"  Bottom layer atoms: {len(bot_super)}")
    print(f"  Top layer atoms: {len(top_super)}")
    
    # Stack layers
    hetero = stack(bot_super, top_super, axis=2, distance=distance)
    hetero.cell[2, 2] = cell_c
    
    print(f"  Total atoms: {len(hetero)}")
    print(f"  Interlayer distance: {distance} Å")
    print(f"  Cell c-axis: {cell_c} Å")
    
    # Analyze composition
    element_counts, element_order = analyze_structure(hetero)
    print(f"  Composition: {dict(element_counts)}")
    print(f"  Element order: {element_order}")
    
    return hetero


def create_single_calculation(
    bot_name,
    top_name,
    bot_file,
    top_file,
    output_base
):
    """
    Create all input files for one combination
    """
    
    print(f"\n{'='*70}")
    print(f"Creating inputs: {bot_name} (bottom) + {top_name} (top)")
    print(f"{'='*70}")
    
    # Create output directory
    dir_name = f"{bot_name}__{top_name}"
    output_dir = os.path.join(output_base, dir_name)
    os.makedirs(output_dir, exist_ok=True)
    
    # Build structure
    hetero = build_heterostructure(
        bot_file,
        top_file,
        supercell,
        interlayer_distance,
        cell_height_c
    )
    
    # Write POSCAR
    poscar_path = os.path.join(output_dir, "POSCAR")
    write(poscar_path, hetero, format="vasp", vasp5=True)
    print(f"  ✓ POSCAR written: {poscar_path}")
    
    # Generate INCAR
    n_cores = pbs_nodes * pbs_ppn
    incar_content = generate_incar(hetero, n_cores, calculation_type)
    incar_path = os.path.join(output_dir, "INCAR")
    with open(incar_path, 'w') as f:
        f.write(incar_content)
    print(f"  ✓ INCAR written: {incar_path}")
    
    # Generate KPOINTS
    kpoints_content, k_mesh = generate_kpoints(
        supercell[0],
        k_mesh_density,
        'relax'
    )
    kpoints_path = os.path.join(output_dir, "KPOINTS")
    with open(kpoints_path, 'w') as f:
        f.write(kpoints_content)
    print(f"  ✓ KPOINTS written: {kpoints_path} ({k_mesh}×{k_mesh}×1)")
    
    # Copy POTCAR
    if os.path.exists(potcar_path):
        os.system(f"cp {potcar_path} {output_dir}/POTCAR")
        print(f"  ✓ POTCAR copied")
    else:
        print(f"  ⚠ Warning: POTCAR not found at {potcar_path}")
    
    # Generate PBS script
    job_name = f"relax_{bot_name}_{top_name}"[:15]  # PBS name limit
    pbs_content = generate_pbs_script(
        job_name=job_name,
        account=pbs_account,
        queue=pbs_queue,
        nodes=pbs_nodes,
        ppn=pbs_ppn,
        gen=pbs_gen,
        walltime=walltime,
        vasp_module=vasp_module,
        email=email,
        notify_mode=notify_mode,
        work_dir=output_dir
    )
    pbs_path = os.path.join(output_dir, "submit.pbs")
    with open(pbs_path, 'w') as f:
        f.write(pbs_content)
    print(f"  ✓ PBS script written: {pbs_path}")
    
    # Write summary
    summary_path = os.path.join(output_dir, "setup_summary.txt")
    with open(summary_path, 'w') as f:
        f.write(f"{'='*70}\n")
        f.write(f"SETUP SUMMARY\n")
        f.write(f"{'='*70}\n\n")
        f.write(f"Bottom layer: {bot_name}\n")
        f.write(f"Top layer: {top_name}\n")
        f.write(f"Supercell: {supercell[0]}×{supercell[1]}×{supercell[2]}\n")
        f.write(f"Total atoms: {len(hetero)}\n")
        f.write(f"Interlayer distance: {interlayer_distance} Å\n")
        f.write(f"Cell c-axis: {cell_height_c} Å\n\n")
        
        element_counts, element_order = analyze_structure(hetero)
        f.write(f"Composition: {dict(element_counts)}\n")
        f.write(f"Element order: {element_order}\n\n")
        
        f.write(f"{'='*70}\n")
        f.write(f"CALCULATION SETTINGS\n")
        f.write(f"{'='*70}\n\n")
        f.write(f"Type: {calculation_type}\n")
        f.write(f"K-mesh: {k_mesh}×{k_mesh}×1 (Gamma-centered)\n")
        f.write(f"Cores: {n_cores} ({pbs_nodes} nodes × {pbs_ppn} ppn)\n")
        f.write(f"Walltime: {walltime}\n")
        f.write(f"Account: {pbs_account}\n")
        f.write(f"Queue: {pbs_queue}\n\n")
        
        f.write(f"{'='*70}\n")
        f.write(f"TO SUBMIT:\n")
        f.write(f"{'='*70}\n\n")
        f.write(f"cd {output_dir}\n")
        f.write(f"qsub submit.pbs\n")
    
    print(f"  ✓ Summary written: {summary_path}")
    print(f"\n  Directory created: {output_dir}")
    
    return output_dir


def main():
    """
    Main execution
    """
    
    print("\n" + "="*70)
    print("JANUS TMD HETEROSTRUCTURE INPUT GENERATOR v2.0")
    print("="*70)
    print("\nImprovements:")
    print("  ✓ Proper MAGMOM for TMDs (Mo/W: 0.5, S/Se: 0.0)")
    print("  ✓ Gamma-centered KPOINTS (mandatory for 2D)")
    print("  ✓ Optimized k-mesh for supercell (3×3×1 or 4×4×1)")
    print("  ✓ ISIF=3 (full cell + ionic relaxation)")
    print("  ✓ Automatic NPAR calculation")
    print("  ✓ vdW-D3(BJ) correction (IVDW=12)")
    print("="*70)
    
    # Check input files
    print("\nChecking input files...")
    missing_files = []
    for name, filepath in structures.items():
        if os.path.exists(filepath):
            print(f"  ✓ {name}: {filepath}")
        else:
            print(f"  ✗ {name}: {filepath} NOT FOUND")
            missing_files.append(filepath)
    
    if missing_files:
        print(f"\n⚠ Error: Missing input files. Please check paths.")
        return
    
    # Create output base directory
    os.makedirs(output_base, exist_ok=True)
    print(f"\nOutput directory: {output_base}")
    
    # Process each combination
    print(f"\n{'='*70}")
    print(f"Generating inputs for {len(combinations_to_run)} combinations")
    print(f"{'='*70}")
    
    created_dirs = []
    for bot_name, top_name in combinations_to_run:
        try:
            output_dir = create_single_calculation(
                bot_name,
                top_name,
                structures[bot_name],
                structures[top_name],
                output_base
            )
            created_dirs.append(output_dir)
        except Exception as e:
            print(f"\n⚠ Error creating {bot_name}+{top_name}: {e}")
            continue
    
    # Final summary
    print("\n" + "="*70)
    print("GENERATION COMPLETE")
    print("="*70)
    print(f"\nCreated {len(created_dirs)} calculation directories:")
    for d in created_dirs:
        print(f"  - {d}")
    
    print(f"\n{'='*70}")
    print("NEXT STEPS:")
    print("="*70)
    print("""
1. Verify POTCAR order matches POSCAR:
   grep TITEL POTCAR
   
2. Submit jobs:
   cd relax_0deg_structures/[directory]
   qsub submit.pbs
   
3. Monitor jobs:
   qstat -u $USER
   
4. Check convergence after completion:
   grep "reached required accuracy" */OUTCAR
   
5. Extract optimized parameters:
   # Use extraction script (to be provided)
    """)
    print("="*70 + "\n")


if __name__ == "__main__":
    main()
