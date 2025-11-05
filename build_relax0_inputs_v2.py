#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
build_relax0_inputs_v2.py
----------------------------------------------------
Generate 0-degree Janus heterostructure (MoSSe–WSSe)
and prepare VASP relaxation input files.

FIXED: Proper POSCAR format with grouped elements and Direct coordinates
----------------------------------------------------
"""

import os
import numpy as np
from ase.io import read
from ase.build import stack
from ase import Atoms
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
# Expected POTCAR order: Mo_sv S Se W_sv S Se (check output)

# Output base directory
output_base = "./relax_0deg_structures"

# PBS and VASP parameters
pbs_account = "cnm84150"
pbs_queue = "batch"
pbs_nodes = 3
pbs_ppn = 16
pbs_gen = "gen6"
walltime = "72:00:00"
vasp_module = "vasp6"
email = "jbaek27@uic.edu"
notify_mode = "abe"

# Calculation settings
calculation_type = "relax"
k_mesh_density = "normal"

# ============================================================
# COMBINATIONS TO CALCULATE
# ============================================================

combinations_to_run = [
    ('MoSSe_S_up', 'WSSe_S_up'),
    ('MoSSe_Se_up', 'WSSe_Se_up'),
    # Add more as needed
]

# ============================================================
# HELPER FUNCTIONS
# ============================================================

def write_poscar_grouped(filename, atoms, comment="Structure"):
    """
    Write POSCAR with properly grouped elements (Direct coordinates)
    Preserves layer structure: bottom layer elements, then top layer elements
    Each layer's elements are grouped separately
    """
    
    # Get all symbols and positions
    symbols = atoms.get_chemical_symbols()
    positions = atoms.get_scaled_positions()
    cell = atoms.get_cell()
    
    # Split into bottom and top layers based on Cartesian z-coordinate
    # More robust for heterostructures with large vacuum
    n_atoms = len(symbols)
    
    # Get Cartesian z-coordinates for splitting
    positions_cart = atoms.get_positions()
    z_coords_cart = positions_cart[:, 2]
    
    # Find the gap between layers in Cartesian coordinates
    z_sorted_indices = np.argsort(z_coords_cart)
    z_sorted_cart = z_coords_cart[z_sorted_indices]
    z_diffs_cart = np.diff(z_sorted_cart)
    max_gap_idx = np.argmax(z_diffs_cart)
    
    # Split point: middle of the largest gap (in Cartesian)
    z_split_cart = (z_sorted_cart[max_gap_idx] + z_sorted_cart[max_gap_idx + 1]) / 2
    
    # Separate atoms into bottom and top layers (Cartesian-based)
    bottom_indices = [i for i, z in enumerate(z_coords_cart) if z < z_split_cart]
    top_indices = [i for i, z in enumerate(z_coords_cart) if z >= z_split_cart]
    
    bottom_symbols = [symbols[i] for i in bottom_indices]
    top_symbols = [symbols[i] for i in top_indices]
    bottom_positions = positions[bottom_indices]  # Keep fractional for POSCAR
    top_positions = positions[top_indices]
    
    print(f"  Layer split at z = {z_split_cart:.4f} Å (Cartesian)")
    print(f"  Split: {len(bottom_symbols)} bottom + {len(top_symbols)} top atoms")
    
    # Debug: show element distribution in each layer
    from collections import Counter
    bottom_counts = Counter(bottom_symbols)
    top_counts = Counter(top_symbols)
    print(f"  Bottom layer: {dict(bottom_counts)}")
    print(f"  Top layer: {dict(top_counts)}")
    
    # Get unique elements in each layer (preserving order)
    def get_unique_ordered(symbol_list):
        unique = []
        seen = set()
        for sym in symbol_list:
            if sym not in seen:
                unique.append(sym)
                seen.add(sym)
        return unique
    
    bottom_elements = get_unique_ordered(bottom_symbols)
    top_elements = get_unique_ordered(top_symbols)
    
    # Build element order: bottom elements + top elements
    element_order = bottom_elements + top_elements
    
    # Group and sort atoms by element within each layer
    # Keep track of which elements came from which layer
    sorted_symbols = []
    sorted_positions = []
    element_counts_ordered = []
    
    # Process bottom layer
    for element in bottom_elements:
        count = 0
        for i, sym in enumerate(bottom_symbols):
            if sym == element:
                sorted_symbols.append(sym)
                sorted_positions.append(bottom_positions[i])
                count += 1
        element_counts_ordered.append(count)
    
    # Process top layer
    for element in top_elements:
        count = 0
        for i, sym in enumerate(top_symbols):
            if sym == element:
                sorted_symbols.append(sym)
                sorted_positions.append(top_positions[i])
                count += 1
        element_counts_ordered.append(count)
    
    sorted_positions = np.array(sorted_positions)
    
    # Write POSCAR
    with open(filename, 'w') as f:
        # Comment line
        f.write(f"{comment}\n")
        
        # Scaling factor
        f.write("1.0\n")
        
        # Lattice vectors
        for vec in cell:
            f.write(f"  {vec[0]:20.16f}  {vec[1]:20.16f}  {vec[2]:20.16f}\n")
        
        # Element names
        f.write(" ".join(element_order) + "\n")
        
        # Element counts
        f.write(" ".join(map(str, element_counts_ordered)) + "\n")
        
        # Direct coordinates
        f.write("Direct\n")
        
        # Write positions with element labels
        atom_index = 0
        for i, elem in enumerate(element_order):
            elem_count = element_counts_ordered[i]
            for j in range(elem_count):
                pos = sorted_positions[atom_index]
                f.write(f"  {pos[0]:20.16f}  {pos[1]:20.16f}  {pos[2]:20.16f}  {elem}\n")
                atom_index += 1
    
    return element_order


def calculate_optimal_kpoints(supercell_size, density='normal'):
    """Calculate optimal k-point mesh for supercell"""
    
    k_base_dict = {
        'coarse': 12,
        'normal': 15,
        'fine': 18
    }
    
    k_base = k_base_dict.get(density, 15)
    k_mesh = max(2, k_base // supercell_size)
    
    return k_mesh


def calculate_npar(n_cores):
    """Calculate optimal NPAR for given number of cores"""
    
    candidates = [4, 6, 8, 12, 16]
    
    for npar in candidates:
        if n_cores % npar == 0:
            return npar
    
    for npar in range(8, 2, -1):
        if n_cores % npar == 0:
            return npar
    
    return 4


def generate_magmom(element_order, element_counts):
    """
    Generate MAGMOM string based on element types and order
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

def generate_incar(element_order, element_counts, n_cores, calculation_type='relax'):
    """Generate INCAR file with proper settings"""
    
    magmom_string = generate_magmom(element_order, element_counts)
    npar = calculate_npar(n_cores)
    
    if calculation_type == 'test':
        nsw = 50
        ediffg = -0.02
    else:
        nsw = 300
        ediffg = -0.015
    
    incar = f"""SYSTEM  = Janus TMD Heterostructure 0deg relax

# === Electronic Structure ===
ENCUT   = 550
PREC    = Accurate
EDIFF   = 1E-6
ALGO    = Normal
NELM    = 150

# === Ionic Relaxation ===
IBRION  = 2
NSW     = {nsw}
EDIFFG  = {ediffg}
POTIM   = 0.5

# === Cell Relaxation ===
ISIF    = 3

# === Spin Polarization ===
ISPIN   = 2
MAGMOM  = {magmom_string}

# === vdW Correction ===
IVDW    = 12

# === Accuracy ===
LASPH   = .TRUE.
LMAXMIX = 4
ADDGRID = .TRUE.

# === k-space Integration ===
ISMEAR  = 0
SIGMA   = 0.05

# === Parallelization ===
NPAR    = {npar}
LREAL   = Auto

# === Output ===
LORBIT  = 11
LWAVE   = .FALSE.
LCHARG  = .FALSE.
NEDOS   = 2000
"""
    
    return incar


def generate_kpoints(supercell_size, density='normal'):
    """Generate KPOINTS file (Gamma-centered)"""
    
    k_mesh = calculate_optimal_kpoints(supercell_size, density)
    
    kpoints = f"""Automatic mesh - RELAX
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
    notify_mode
):
    """Generate PBS job script"""
    
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

echo "========================================"
echo "Job started on $(hostname) at $(date)"
echo "Job ID: $PBS_JOBID"
echo "========================================"

module load {vasp_module}
module list

cd $PBS_O_WORKDIR

echo "Running in: $(pwd)"
cat $PBS_NODEFILE

echo "========================================"
echo "INCAR settings:"
grep -E "ISIF|ISPIN|IVDW|ENCUT|NSW" INCAR
echo ""
echo "KPOINTS:"
head -4 KPOINTS
echo "========================================"

echo "Starting VASP at $(date)"
time mpirun -np {n_cores} vasp_std > vasp.log 2>&1

if [ -f CONTCAR ]; then
    echo "Calculation completed at $(date)"
    if grep -q "reached required accuracy" OUTCAR; then
        echo "✓ Convergence achieved!"
    else
        echo "⚠ Check convergence in OUTCAR"
    fi
    
    if [ -f OSZICAR ]; then
        echo ""
        echo "Final energies:"
        tail -5 OSZICAR
    fi
else
    echo "⚠ CONTCAR not found - calculation may have failed"
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
    """Build heterostructure from two layer files - MANUAL STACKING"""
    
    print(f"  Reading bottom: {os.path.basename(bot_file)}")
    print(f"  Reading top: {os.path.basename(top_file)}")
    
    bot = read(bot_file)
    top = read(top_file)
    
    # Build supercell
    bot_super = bot.repeat(supercell)
    top_super = top.repeat(supercell)
    
    print(f"  Supercell: {supercell[0]}×{supercell[1]}×{supercell[2]}")
    print(f"  Bottom atoms: {len(bot_super)}")
    print(f"  Top atoms: {len(top_super)}")
    
    # MANUAL STACKING with proper periodic boundary handling
    # Get fractional coordinates first to handle periodicity
    bot_frac = bot_super.get_scaled_positions()
    top_frac = top_super.get_scaled_positions()
    
    # Wrap z-coordinates to [0, 1) and then to compact representation
    # For Janus TMD: atoms with z > 0.5 are actually below (periodic)
    bot_z_frac = bot_frac[:, 2].copy()
    top_z_frac = top_frac[:, 2].copy()
    
    # Wrap to [-0.5, 0.5) to identify atoms that should be below
    bot_z_frac = bot_z_frac % 1.0
    bot_z_frac[bot_z_frac > 0.5] -= 1.0
    
    top_z_frac = top_z_frac % 1.0
    top_z_frac[top_z_frac > 0.5] -= 1.0
    
    # Convert to Cartesian (Å) using original cell
    bot_cell_c = bot_super.cell[2, 2]
    top_cell_c = top_super.cell[2, 2]
    
    bot_z_cart = bot_z_frac * bot_cell_c
    top_z_cart = top_z_frac * top_cell_c
    
    print(f"  DEBUG: After wrapping:")
    print(f"    Bot z: {bot_z_cart.min():.4f} to {bot_z_cart.max():.4f} Å")
    print(f"    Top z: {top_z_cart.min():.4f} to {top_z_cart.max():.4f} Å")
    
    # Shift both to start at z=0
    bot_z_cart -= bot_z_cart.min()
    bot_z_max = bot_z_cart.max()
    bot_thickness = bot_z_max
    
    top_z_cart -= top_z_cart.min()
    top_z_max_rel = top_z_cart.max()
    top_thickness = top_z_max_rel
    
    print(f"  Bottom layer: 0.000 to {bot_z_max:.4f} Å (thickness: {bot_thickness:.4f} Å)")
    
    # Position top layer above bottom
    top_z_start = bot_z_max + distance
    top_z_cart += top_z_start
    top_z_max = top_z_cart.max()
    
    print(f"  Top layer: {top_z_start:.4f} to {top_z_max:.4f} Å (thickness: {top_thickness:.4f} Å)")
    print(f"  Interlayer gap: {distance:.2f} Å")
    
    # Get full Cartesian positions and update z
    bot_positions = bot_super.get_positions()
    top_positions = top_super.get_positions()
    
    bot_positions[:, 2] = bot_z_cart
    top_positions[:, 2] = top_z_cart
    
    # Combine
    all_symbols = bot_super.get_chemical_symbols() + top_super.get_chemical_symbols()
    all_positions = np.vstack([bot_positions, top_positions])
    
    # Create new cell
    cell = bot_super.get_cell().copy()
    cell[2, 2] = cell_c
    
    print(f"  New cell c-axis: {cell_c:.2f} Å")
    
    # Create heterostructure
    hetero = Atoms(
        symbols=all_symbols,
        positions=all_positions,
        cell=cell,
        pbc=[True, True, True]
    )
    
    # Verify
    frac_pos = hetero.get_scaled_positions()
    cart_pos = hetero.get_positions()
    
    print(f"  Total atoms: {len(hetero)}")
    print(f"  Final z range: {cart_pos[:, 2].min():.4f} to {cart_pos[:, 2].max():.4f} Å")
    print(f"  Final z range (fractional): {frac_pos[:, 2].min():.6f} to {frac_pos[:, 2].max():.6f}")
    
    return hetero


def create_single_calculation(
    bot_name,
    top_name,
    bot_file,
    top_file,
    output_base
):
    """Create all input files for one combination"""
    
    print(f"\n{'='*70}")
    print(f"Creating: {bot_name} + {top_name}")
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
    
    # Write POSCAR with grouped elements
    poscar_path = os.path.join(output_dir, "POSCAR")
    comment = f"{bot_name} + {top_name} - 5x5x1 supercell"
    element_order = write_poscar_grouped(poscar_path, hetero, comment)
    print(f"  ✓ POSCAR written: {poscar_path}")
    print(f"    Element order: {element_order}")
    
    # Get element counts
    symbols = hetero.get_chemical_symbols()
    element_counts = Counter(symbols)
    print(f"    Element counts: {dict(element_counts)}")
    
    # Generate INCAR
    n_cores = pbs_nodes * pbs_ppn
    incar_content = generate_incar(element_order, element_counts, n_cores, calculation_type)
    incar_path = os.path.join(output_dir, "INCAR")
    with open(incar_path, 'w') as f:
        f.write(incar_content)
    print(f"  ✓ INCAR written")
    
    # Generate KPOINTS
    kpoints_content, k_mesh = generate_kpoints(supercell[0], k_mesh_density)
    kpoints_path = os.path.join(output_dir, "KPOINTS")
    with open(kpoints_path, 'w') as f:
        f.write(kpoints_content)
    print(f"  ✓ KPOINTS written ({k_mesh}×{k_mesh}×1)")
    
    # Copy POTCAR if exists
    if os.path.exists(potcar_path):
        os.system(f"cp {potcar_path} {output_dir}/POTCAR")
        print(f"  ✓ POTCAR copied")
    else:
        print(f"  ⚠ POTCAR not found (will generate on Carbon)")
    
    # Generate PBS script
    job_name = f"rx_{bot_name[:4]}_{top_name[:4]}"[:15]
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
        notify_mode=notify_mode
    )
    pbs_path = os.path.join(output_dir, "submit.pbs")
    with open(pbs_path, 'w') as f:
        f.write(pbs_content)
    print(f"  ✓ PBS script written")
    
    # Write summary
    summary_path = os.path.join(output_dir, "setup_summary.txt")
    with open(summary_path, 'w') as f:
        f.write(f"{'='*70}\n")
        f.write(f"SETUP SUMMARY\n")
        f.write(f"{'='*70}\n\n")
        f.write(f"Bottom: {bot_name}\n")
        f.write(f"Top: {top_name}\n")
        f.write(f"Supercell: {supercell[0]}×{supercell[1]}×{supercell[2]}\n")
        f.write(f"Total atoms: {len(hetero)}\n")
        f.write(f"Interlayer distance: {interlayer_distance} Å (initial)\n")
        f.write(f"Cell c-axis: {cell_height_c} Å\n\n")
        
        f.write(f"Element order: {' '.join(element_order)}\n")
        f.write(f"Element counts: {' '.join(str(element_counts[e]) for e in element_order)}\n\n")
        
        f.write(f"POTCAR order needed:\n")
        for i, elem in enumerate(element_order, 1):
            if elem in ['Mo', 'W']:
                f.write(f"  {i}. {elem}_sv\n")
            else:
                f.write(f"  {i}. {elem}\n")
        f.write(f"\n{'='*70}\n")
        f.write(f"TO SUBMIT:\n")
        f.write(f"cd {output_dir}\n")
        f.write(f"qsub submit.pbs\n")
    
    print(f"  ✓ Summary written")
    print(f"\n  Directory: {output_dir}")
    
    return output_dir


def main():
    """Main execution"""
    
    print("\n" + "="*70)
    print("JANUS TMD HETEROSTRUCTURE INPUT GENERATOR v2.1")
    print("="*70)
    print("\nFEATURES:")
    print("  ✓ Proper POSCAR format (Direct coordinates)")
    print("  ✓ Grouped elements (matching reference)")
    print("  ✓ Proper MAGMOM (Mo/W: 0.5, S/Se: 0.0)")
    print("  ✓ Gamma-centered KPOINTS")
    print("  ✓ ISIF=3 full relaxation")
    print("  ✓ vdW-D3(BJ) correction")
    print("="*70)
    
    # Check input files
    print("\nChecking input files...")
    missing_files = []
    for name, filepath in structures.items():
        if os.path.exists(filepath):
            print(f"  ✓ {name}")
        else:
            print(f"  ✗ {name}: NOT FOUND")
            missing_files.append(filepath)
    
    if missing_files:
        print(f"\n⚠ Error: Missing files")
        return
    
    # Create output directory
    os.makedirs(output_base, exist_ok=True)
    print(f"\nOutput: {output_base}")
    
    # Process combinations
    print(f"\n{'='*70}")
    print(f"Generating {len(combinations_to_run)} combinations")
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
            print(f"\n⚠ Error: {e}")
            import traceback
            traceback.print_exc()
            continue
    
    # Final summary
    print("\n" + "="*70)
    print("GENERATION COMPLETE")
    print("="*70)
    print(f"\nCreated {len(created_dirs)} directories")
    
    print(f"\n{'='*70}")
    print("NEXT STEPS:")
    print("="*70)
    print("""
1. Verify POSCAR format:
   head -10 relax_0deg_structures/*/POSCAR
   
2. Transfer to Carbon:
   rsync -avz relax_0deg_structures/ user@carbon:~/project/
   
3. Generate POTCAR on Carbon:
   ./generate_potcar.sh
   
4. Verify POTCAR order:
   grep TITEL */POTCAR | head -6
   
5. Submit jobs:
   ./submit_all_jobs.sh
    """)
    print("="*70 + "\n")


if __name__ == "__main__":
    main()
