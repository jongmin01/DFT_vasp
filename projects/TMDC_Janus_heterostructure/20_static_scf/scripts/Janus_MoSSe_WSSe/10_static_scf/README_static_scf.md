# Static SCF Calculation Guide

## Overview

Static SCF (Self-Consistent Field) calculation performs precise electronic structure calculation on **fixed, optimized structures**.

**Purpose:**
- Accurate total energy
- Band structure preparation
- Density of States (DOS)
- Charge density analysis
- Optical properties

**Key difference from relaxation:**
- NSW = 0 (no ionic steps)
- IBRION = -1 (no relaxation)
- Higher accuracy (EDIFF = 1E-06)

## Files Generated

### Input Files:
1. **[INCAR_static_scf](computer:///mnt/user-data/outputs/INCAR_static_scf)** - Static calculation settings
2. **[KPOINTS_static_scf](computer:///mnt/user-data/outputs/KPOINTS_static_scf)** - K-point mesh (7×7×1)
3. **[run_static_scf.pbs](computer:///mnt/user-data/outputs/run_static_scf.pbs)** - PBS job script
4. **[setup_static_scf.sh](computer:///mnt/user-data/outputs/setup_static_scf.sh)** - Automated setup
5. **[analyze_static_scf.sh](computer:///mnt/user-data/outputs/analyze_static_scf.sh)** - Results analysis

## Quick Start

### Method 1: Automated Setup (Recommended)

```bash
# 1. Download all files to your calculation directory
cd /path/to/calculations

# 2. Make scripts executable
chmod +x setup_static_scf.sh analyze_static_scf.sh

# 3. Run setup (will ask about job submission)
./setup_static_scf.sh
```

This will:
- Create `static_scf/` subdirectories
- Copy relaxed structures (CONTCAR_final → POSCAR)
- Set up all input files
- Optionally submit jobs

### Method 2: Manual Setup

```bash
# For each directory
cd MoSSe_Se_up__WSSe_S_up

# Create subdirectory
mkdir static_scf
cd static_scf

# Copy relaxed structure
cp ../CONTCAR_final POSCAR

# Copy input files
cp ../../INCAR_static_scf INCAR
cp ../../KPOINTS_static_scf KPOINTS
cp ../../run_static_scf.pbs .
cp ../POTCAR .

# Submit
qsub run_static_scf.pbs
```

## Directory Structure

After setup:
```
MoSSe_Se_up__WSSe_S_up/
├── stage1_results/          # Relaxation results
│   └── CONTCAR
├── CONTCAR_final            # Best relaxed structure
├── POSCAR_for_static_scf    # Ready for static
└── static_scf/              # NEW - Static calculation
    ├── POSCAR               # Relaxed structure
    ├── INCAR                # Static settings
    ├── KPOINTS              # 7×7×1 mesh
    ├── POTCAR               # Pseudopotentials
    └── run_static_scf.pbs   # Job script
```

## INCAR Settings Explained

### Key Parameters:

```bash
# No relaxation
IBRION = -1          # No ionic movement
NSW = 0              # Zero ionic steps

# High accuracy
EDIFF = 1E-06        # Tight electronic convergence

# Output control
LORBIT = 11          # Write DOSCAR with projections
NEDOS = 3000         # DOS energy resolution
LWAVE = .TRUE.       # Write WAVECAR (for bands)
LCHARG = .TRUE.      # Write CHGCAR
LAECHG = .TRUE.      # Write core charges
```

### Compared to Relaxation:

| Parameter | Relaxation | Static SCF |
|-----------|-----------|------------|
| IBRION | 2 | -1 |
| NSW | 150-300 | 0 |
| EDIFF | 5E-06 | 1E-06 |
| LORBIT | - | 11 |
| LWAVE | .FALSE. | .TRUE. |

## Expected Resources

### Computation Time:
- **Per structure**: 8-12 hours
- **Total (4 structures)**: ~12 hours (parallel)
- **Nodes**: 2 nodes, 32 cores
- **Walltime**: 24 hours (safe)

### Disk Space:
```
OUTCAR:   ~50 MB
WAVECAR:  ~2-5 GB  ← Large!
CHGCAR:   ~1-2 GB
DOSCAR:   ~5 MB
EIGENVAL: ~20 MB
Total:    ~4-8 GB per structure
```

**Important**: WAVECAR is large - clean up after use if space limited!

## Monitoring

### Real-time monitoring:
```bash
# Job status
qstat -u $USER

# Watch output
tail -f MoSSe_Se_up__WSSe_S_up/static_scf/*.o*

# Electronic iterations
grep "DAV:" MoSSe_Se_up__WSSe_S_up/static_scf/OSZICAR | tail -10
```

### Progress check:
```bash
# Quick status
for dir in MoSSe_*/static_scf; do
    echo "$dir: $(grep -c 'DAV:' $dir/OSZICAR 2>/dev/null || echo 0) iterations"
done
```

### Analyze results:
```bash
./analyze_static_scf.sh
```

## Expected Output

### Convergence:
```
Electronic iterations: 20-40 (typical)
Convergence: "reached required accuracy"
```

### Energy:
```
Similar to relaxation final energy
Difference: < 0.001 eV (should be very small)
```

### Files created:
- ✅ OUTCAR - Main output
- ✅ DOSCAR - Density of states
- ✅ EIGENVAL - Eigenvalues
- ✅ CHGCAR - Charge density
- ✅ WAVECAR - Wavefunctions
- ✅ vasprun.xml - XML output

## Next Steps After Static SCF

### 1. DOS Analysis

```bash
# Using p4vasp (if available)
p4v MoSSe_Se_up__WSSe_S_up/static_scf/vasprun.xml

# Or extract manually
cd MoSSe_Se_up__WSSe_S_up/static_scf
# Use vaspkit, pymatgen, or custom scripts
```

### 2. Band Structure Calculation

Requires additional calculation along high-symmetry k-path:
- Will need separate KPOINTS with path
- Use ICHARG = 11 (read charge density)

### 3. Charge Analysis

```bash
# Bader analysis
bader CHGCAR -ref CHGCAR_sum

# Or use VASP tools
```

### 4. Optical Properties (Optional)

Add to INCAR:
```bash
LOPTICS = .TRUE.
NBANDS = 400  # More bands
```

## Troubleshooting

### Job fails immediately:
```bash
# Check files exist
ls -lh MoSSe_*/static_scf/POSCAR

# Verify POTCAR
grep TITEL MoSSe_Se_up__WSSe_S_up/static_scf/POTCAR
```

### Doesn't converge:
```bash
# Increase NELM
NELM = 200

# Try different algorithm
ALGO = Fast  # or All
```

### WAVECAR too large:
```bash
# After calculation, can delete if not needed
rm MoSSe_*/static_scf/WAVECAR

# Or compress
gzip MoSSe_*/static_scf/WAVECAR
```

### Out of disk space:
```bash
# Check usage
du -sh MoSSe_*/static_scf

# Remove unnecessary files
rm MoSSe_*/static_scf/CHG  # Can regenerate from CHGCAR
```

## Validation Checklist

After calculation completes:

- [ ] All 4 structures converged?
- [ ] Energy similar to relaxation (< 0.01 eV difference)?
- [ ] DOSCAR file exists and readable?
- [ ] EIGENVAL file exists?
- [ ] No errors in OUTCAR?
- [ ] Fermi level reasonable?

## Summary

**What you get:**
- ✅ Precise electronic structure
- ✅ DOS for analysis
- ✅ Wavefunctions for bands
- ✅ Charge density
- ✅ Ready for advanced analysis

**Time investment:**
- Setup: 5 minutes
- Calculation: 12 hours
- Analysis: Variable

**Next analysis options:**
1. DOS plotting
2. Band structure (needs additional calculation)
3. Charge analysis (Bader)
4. Comparison between stackings

---

**Ready to start? Run `./setup_static_scf.sh` now!**
