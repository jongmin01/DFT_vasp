#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
generate_pbs_scripts.py
Generate PBS job scripts for all heterostructure directories
Python 2/3 compatible
"""

from __future__ import print_function
import os
import sys

# ============================================================
# PBS SETTINGS - MODIFY THESE FOR YOUR SYSTEM
# ============================================================

# Account and queue settings
pbs_account = "cnm84150"
pbs_queue = "batch"

# Node settings
pbs_nodes = 3
pbs_ppn = 16
pbs_gen = "gen6"  # gen6 has 192GB RAM (needed for 150 atoms)

# Walltime
walltime = "72:00:00"

# VASP module
vasp_module = "vasp5"  # or "vasp6" or specific version like "vasp/5.4.4"

# Email notifications
email = "jbaek27@uic.edu"
notify_mode = "abe"  # a=abort, b=begin, e=end

# Directories to process (auto-detect or specify)
auto_detect = True  # Set to False to use manual list below
manual_dirs = [
    'MoSSe_S_up__WSSe_S_up',
    'MoSSe_S_up__WSSe_Se_up',
    'MoSSe_Se_up__WSSe_S_up',
    'MoSSe_Se_up__WSSe_Se_up',
]

# ============================================================
# PBS SCRIPT TEMPLATE
# ============================================================

def generate_pbs_script(job_name, n_cores):
    """Generate PBS job script content"""
    
    script = """#!/bin/bash
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
# Generated automatically for Janus heterostructure
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

# Check required files
echo "========================================"
echo "Checking input files..."
for file in POSCAR INCAR KPOINTS POTCAR; do
    if [ -f $file ]; then
        echo "  OK: $file ($(stat -f%z $file 2>/dev/null || stat -c%s $file) bytes)"
    else
        echo "  ERROR: $file not found!"
        exit 1
    fi
done
echo "========================================"

# Print INCAR settings
echo ""
echo "INCAR settings:"
grep -E "PREC|ENCUT|ALGO|IBRION|ISIF|NSW|EDIFFG" INCAR | head -10
echo ""

# Print KPOINTS
echo "KPOINTS:"
head -4 KPOINTS
echo ""

# Print structure info
echo "Structure info:"
head -7 POSCAR | tail -2
echo "========================================"
echo ""

# Run VASP
echo "Starting VASP calculation at $(date)"
echo "Command: mpirun -np {n_cores} vasp_std"
echo ""

time mpirun -np {n_cores} vasp_std > vasp.log 2>&1

# Check completion
echo ""
echo "========================================"
echo "Job finished at $(date)"
echo "========================================"
echo ""

if [ -f CONTCAR ]; then
    echo "CONTCAR found - calculation completed"
    echo ""
    
    # Check convergence
    if grep -q "reached required accuracy" OUTCAR; then
        echo "*** CONVERGENCE ACHIEVED ***"
    else
        echo "WARNING: Check convergence in OUTCAR"
    fi
    echo ""
    
    # Extract final energy and forces
    if [ -f OSZICAR ]; then
        echo "Final ionic steps (from OSZICAR):"
        tail -5 OSZICAR
        echo ""
    fi
    
    # Check for errors
    if grep -q "ERROR" OUTCAR 2>/dev/null; then
        echo "WARNING: Errors found in OUTCAR:"
        grep "ERROR" OUTCAR | head -5
        echo ""
    fi
    
    # Summary
    echo "Job summary:"
    echo "  Working directory: $PBS_O_WORKDIR"
    echo "  Job ID: $PBS_JOBID"
    echo "  Walltime used: $(qstat -f $PBS_JOBID 2>/dev/null | grep resources_used.walltime || echo 'N/A')"
    
else
    echo "*** WARNING: CONTCAR not found ***"
    echo "Calculation may have failed or is incomplete"
    echo ""
    echo "Check these files for errors:"
    echo "  - OUTCAR (last 50 lines)"
    echo "  - vasp.log"
    echo "  - Job output file"
    echo ""
    
    if [ -f OUTCAR ]; then
        echo "Last 20 lines of OUTCAR:"
        tail -20 OUTCAR
    fi
fi

echo "========================================"
echo "End of job script"
echo "========================================"
""".format(
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
        n_cores=n_cores
    )
    
    return script

# ============================================================
# MAIN FUNCTION
# ============================================================

def main():
    """Generate PBS scripts for all directories"""
    
    n_cores = pbs_nodes * pbs_ppn
    
    print("=" * 60)
    print("PBS Job Script Generator")
    print("=" * 60)
    print("Settings:")
    print("  Account: " + pbs_account)
    print("  Queue: " + pbs_queue)
    print("  Nodes: {0} x {1} cores = {2} cores total".format(
        pbs_nodes, pbs_ppn, n_cores))
    print("  Generation: " + pbs_gen)
    print("  Walltime: " + walltime)
    print("  VASP module: " + vasp_module)
    print("  Email: " + email)
    print()
    
    # Get directories to process
    if auto_detect:
        # Find all MoSSe_*__WSSe_* directories
        all_dirs = [d for d in os.listdir('.') if os.path.isdir(d)]
        target_dirs = [d for d in all_dirs if d.startswith('MoSSe_') and '__WSSe_' in d]
        
        if not target_dirs:
            print("No heterostructure directories found!")
            print("Looking for: MoSSe_*__WSSe_*")
            print()
            print("Available directories:")
            for d in sorted(all_dirs):
                print("  " + d)
            sys.exit(1)
    else:
        target_dirs = manual_dirs
    
    print("Target directories ({0}):".format(len(target_dirs)))
    for d in sorted(target_dirs):
        print("  " + d)
    print()
    
    # Generate PBS scripts
    print("Generating PBS scripts...")
    print()
    
    success_count = 0
    for dirname in sorted(target_dirs):
        if not os.path.isdir(dirname):
            print("  SKIP: {0} (not found)".format(dirname))
            continue
        
        # Check if required files exist
        required_files = ['POSCAR', 'INCAR', 'KPOINTS']
        missing_files = [f for f in required_files if not os.path.exists(os.path.join(dirname, f))]
        
        if missing_files:
            print("  SKIP: {0} (missing: {1})".format(
                dirname, ', '.join(missing_files)))
            continue
        
        # Generate job name (shorten if too long)
        job_name = dirname
        if len(job_name) > 15:
            # PBS job names should be <= 15 chars
            # Extract key parts: MoSSe_Se_up__WSSe_S_up -> Mo_Se_W_S
            parts = job_name.split('__')
            if len(parts) == 2:
                # MoSSe_Se_up -> Mo_Se, WSSe_S_up -> W_S
                bot = parts[0].split('_')  # ['MoSSe', 'Se', 'up']
                top = parts[1].split('_')  # ['WSSe', 'S', 'up']
                if len(bot) >= 2 and len(top) >= 2:
                    job_name = "{0}_{1}__{2}_{3}".format(
                        bot[0][:2], bot[1], top[0][:1], top[1])
        
        # Generate script
        script_content = generate_pbs_script(job_name, n_cores)
        script_path = os.path.join(dirname, 'run.pbs')
        
        with open(script_path, 'w') as f:
            f.write(script_content)
        
        # Make executable
        os.chmod(script_path, 0o755)
        
        print("  OK: {0}/run.pbs".format(dirname))
        success_count += 1
    
    print()
    print("=" * 60)
    print("Generated {0} PBS scripts".format(success_count))
    print("=" * 60)
    print()
    
    if success_count > 0:
        print("Next steps:")
        print()
        print("1. Review a sample script:")
        print("   cat {0}/run.pbs".format(sorted(target_dirs)[0]))
        print()
        print("2. Submit jobs:")
        print("   for d in MoSSe_*__WSSe_*/; do")
        print("       cd \"$d\"")
        print("       qsub run.pbs")
        print("       cd ..")
        print("   done")
        print()
        print("3. Or submit individually:")
        print("   cd MoSSe_Se_up__WSSe_S_up")
        print("   qsub run.pbs")
        print()
        print("4. Check job status:")
        print("   qstat -u $USER")
        print("=" * 60)
    else:
        print("ERROR: No PBS scripts were generated!")
        print("Check that directories exist and contain required input files")

if __name__ == '__main__':
    main()
