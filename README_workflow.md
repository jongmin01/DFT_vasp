# Janus TMD Heterostructure VASP Calculation Workflow

Complete workflow for building and running VASP calculations on Janus heterostructures (MoSSe-WSSe).

## Quick Start

```bash
# 1. Structure generation
python build_heterostructure_from_seeds.py

# 2. POTCAR generation
python generate_potcar_from_Z.py

# 3. Copy POTCAR to all directories
for d in MoSSe_*__WSSe_*/; do cp POTCAR "$d/"; done

# 4. Generate PBS job scripts
python generate_pbs_scripts.py

# 5. Submit all jobs
bash submit_all_jobs.sh

# 6. Monitor progress
bash check_calculations.sh
```

## Detailed Steps

### Step 1: Prepare Seed Files

Make sure you have these files in your working directory:
- `MoSeS_S_up.vasp`
- `MoSeS_Se_up.vasp`
- `WSeS_S_up.vasp`
- `WSeS_Se_up.vasp`

### Step 2: Generate Heterostructures

```bash
python build_heterostructure_from_seeds.py
```

This creates 4 directories:
- `MoSSe_S_up__WSSe_S_up/`
- `MoSSe_S_up__WSSe_Se_up/`
- `MoSSe_Se_up__WSSe_S_up/`
- `MoSSe_Se_up__WSSe_Se_up/`

Each directory contains:
- `POSCAR` (150 atoms, 5x5 supercell)
- `INCAR` (relaxation settings)
- `KPOINTS` (2x2x1 mesh)

### Step 3: Generate POTCAR

#### Option A: Using Python script
```bash
# Edit generate_potcar_from_Z.py to set correct POTCAR_BASE path
# Default: /opt/soft/vasp-pot/ase/potpaw_PBE

python generate_potcar_from_Z.py
```

#### Option B: Manual command
```bash
zcat /opt/soft/vasp-pot/ase/potpaw_PBE/Mo_pv/POTCAR.Z \
     /opt/soft/vasp-pot/ase/potpaw_PBE/W_pv/POTCAR.Z \
     /opt/soft/vasp-pot/ase/potpaw_PBE/S/POTCAR.Z \
     /opt/soft/vasp-pot/ase/potpaw_PBE/Se/POTCAR.Z > POTCAR
```

#### Verify POTCAR
```bash
grep TITEL POTCAR
# Should show 4 lines:
#   TITEL  = PAW_PBE Mo_pv 08Apr2002
#   TITEL  = PAW_PBE W_pv 06Sep2000
#   TITEL  = PAW_PBE S 06Sep2000
#   TITEL  = PAW_PBE Se 06Sep2000
```

### Step 4: Copy POTCAR to All Directories

```bash
for d in MoSSe_*__WSSe_*/; do
    cp POTCAR "$d/"
    echo "Copied to $d"
done
```

### Step 5: Generate PBS Job Scripts

```bash
# Edit generate_pbs_scripts.py to set your settings:
# - pbs_account
# - pbs_nodes, pbs_ppn, pbs_gen
# - email
# - walltime

python generate_pbs_scripts.py
```

This creates `vasp.job` in each directory.

### Step 6: Review Job Script (Optional)

```bash
cat MoSSe_Se_up__WSSe_S_up/run.pbs
```

### Step 7: Submit Jobs

#### Option A: Submit all at once
```bash
chmod +x submit_all_jobs.sh
./submit_all_jobs.sh
```

#### Option B: Submit individually
```bash
cd MoSSe_Se_up__WSSe_S_up
qsub run.pbs
cd ..
```

### Step 8: Monitor Calculations

```bash
# Check overall status
chmod +x check_calculations.sh
./check_calculations.sh

# Check queue
qstat -u $USER

# Follow output of specific job
tail -f MoSSe_Se_up__WSSe_S_up/*.o<jobid>

# Check convergence
grep "reached required accuracy" MoSSe_*/OUTCAR
```

## File Structure

After completing all steps:

```
.
├── MoSeS_S_up.vasp                    # Seed files
├── MoSeS_Se_up.vasp
├── WSeS_S_up.vasp
├── WSeS_Se_up.vasp
├── POTCAR                             # Master POTCAR
├── build_heterostructure_from_seeds.py
├── generate_potcar_from_Z.py
├── generate_pbs_scripts.py
├── submit_all_jobs.sh
├── check_calculations.sh
│
├── MoSSe_S_up__WSSe_S_up/
│   ├── POSCAR
│   ├── INCAR
│   ├── KPOINTS
│   ├── POTCAR
│   ├── run.pbs
│   └── [output files after run]
│
├── MoSSe_S_up__WSSe_Se_up/
│   └── ...
│
├── MoSSe_Se_up__WSSe_S_up/
│   └── ...
│
└── MoSSe_Se_up__WSSe_Se_up/
    └── ...
```

## Settings Summary

### Structure Settings
- Supercell: 5x5x1
- Total atoms: 150 per structure
- Interlayer distance: 6.8 Angstrom (initial guess)
- Cell height: 30 Angstrom

### INCAR Settings
- PREC = Normal
- ENCUT = 400 eV
- ALGO = Normal
- IBRION = 2 (CG relaxation)
- ISIF = 2 (relax ions only)
- NSW = 200
- EDIFFG = -0.02 eV/Angstrom
- ISMEAR = 0, SIGMA = 0.1
- AMIX = 0.2, BMIX = 0.0001 (conservative mixing)

### PBS Settings
- Nodes: 3 x 16 cores = 48 cores
- Generation: gen6 (192GB RAM)
- Walltime: 72:00:00
- Queue: batch

## Troubleshooting

### Issue: POTCAR species mismatch
```bash
# Check species order in POSCAR
head -7 MoSSe_*/POSCAR | grep -A1 "Direct"

# Should be: Mo  W  S  Se
# POTCAR order must match
```

### Issue: Job fails immediately
```bash
# Check OUTCAR for errors
tail -50 <dir>/OUTCAR

# Common errors:
# - POTCAR missing or wrong order
# - Insufficient memory
# - Module not loaded
```

### Issue: Calculation not converging
```bash
# Check OSZICAR
tail -20 <dir>/OSZICAR

# If energy oscillating:
# - Increase AMIX to 0.1
# - Try ALGO = Fast or All
# - Check initial structure
```

### Issue: All atoms have same z-coordinate
```bash
# This was the original bug - should be fixed
# Verify:
awk 'NR>8 && NF==3 {print $3}' <dir>/POSCAR | sort -u

# Should show 5-6 different z values
```

## Expected Runtime

- With 48 cores (3 nodes x 16 cores)
- 150 atoms
- 200 ionic steps max
- Estimated: 24-48 hours per structure

## Analysis After Completion

```bash
# Extract final energies
for d in MoSSe_*__WSSe_*/; do
    echo "$d:"
    grep "free  energy" "$d/OUTCAR" | tail -1
done

# Check convergence
for d in MoSSe_*__WSSe_*/; do
    if grep -q "reached required accuracy" "$d/OUTCAR"; then
        echo "$d: CONVERGED"
    else
        echo "$d: NOT CONVERGED"
    fi
done

# Extract final structures
for d in MoSSe_*__WSSe_*/; do
    cp "$d/CONTCAR" "${d%/}_final.vasp"
done
```

## Support

For issues or questions:
1. Check OUTCAR for error messages
2. Review PBS job output files (.o<jobid>)
3. Verify input files (POSCAR, INCAR, KPOINTS, POTCAR)
4. Check module availability: `module avail vasp`

## References

- VASP Manual: https://www.vasp.at/wiki/
- Carbon Cluster Wiki: https://wiki.anl.gov/cnm/HPC/
- PBS Documentation: Torque Resource Manager
