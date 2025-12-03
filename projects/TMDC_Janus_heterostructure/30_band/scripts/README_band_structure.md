# Band Structure Calculation Guide (SOC OFF)

## Overview

Band structure calculation provides:
- Energy dispersion E(k) along high-symmetry path
- Direct vs indirect band gap determination
- VBM and CBM locations in k-space
- Effective mass (from band curvature)

## Files Generated

1. **[INCAR_band_soc_off](computer:///mnt/user-data/outputs/INCAR_band_soc_off)** - Band calculation settings
2. **[KPOINTS_band](computer:///mnt/user-data/outputs/KPOINTS_band)** - k-path (Γ-M-K-Γ)
3. **[run_band_soc_off.pbs](computer:///mnt/user-data/outputs/run_band_soc_off.pbs)** - PBS script
4. **[setup_band_structure.sh](computer:///mnt/user-data/outputs/setup_band_structure.sh)** - Automated setup
5. **[plot_band_structure.py](computer:///mnt/user-data/outputs/plot_band_structure.py)** - Plotting script

## Quick Start

### Step 1: Setup

```bash
# Upload files to calculation directory
cd /path/to/your/calculations

# Make executable
chmod +x setup_band_structure.sh

# Run setup
./setup_band_structure.sh
```

This will:
- Create `MoSSe_Se_up__WSSe_S_up/band_soc_off/` directory
- Copy POSCAR, POTCAR, CHGCAR from static_scf
- Setup INCAR with ICHARG=11
- Setup KPOINTS with Γ-M-K-Γ path
- Ask if you want to submit job

### Step 2: Monitor

```bash
# Check job status
qstat -u $USER

# Monitor progress
tail -f MoSSe_Se_up__WSSe_S_up/band_soc_off/*.o*

# Check convergence
grep "F=" MoSSe_Se_up__WSSe_S_up/band_soc_off/OSZICAR
```

### Step 3: Analysis (after completion)

```bash
# Plot band structure
cd MoSSe_Se_up__WSSe_S_up/band_soc_off
python ../../plot_band_structure.py EIGENVAL

# Or specify path
python plot_band_structure.py MoSSe_Se_up__WSSe_S_up/band_soc_off/EIGENVAL
```

## Directory Structure

After setup:
```
MoSSe_Se_up__WSSe_S_up/
├── stage1_results/     # Relaxation
├── static_scf/         # DOS, charge density
│   ├── POSCAR
│   ├── POTCAR
│   └── CHGCAR         # ← Used for band
└── band_soc_off/       # NEW - Band structure
    ├── POSCAR         # From static_scf
    ├── INCAR          # ICHARG=11
    ├── KPOINTS        # Γ-M-K-Γ path
    ├── POTCAR         # From static_scf
    ├── CHGCAR         # From static_scf
    └── run_band_soc_off.pbs
```

## INCAR Settings Explained

### Key Parameters:

```bash
# Reuse charge density
ICHARG = 11            # Read CHGCAR, fixed charge density

# No relaxation
IBRION = -1
NSW = 0

# Spin-polarized (same as static SCF)
ISPIN = 2

# vdW included
IVDW = 12

# Output
LORBIT = 11            # For orbital projection
NEDOS = 3000
```

### Differences from Static SCF:

| Parameter | Static SCF | Band Structure |
|-----------|-----------|----------------|
| ICHARG | 2 (charge) | 11 (read CHGCAR) |
| KPOINTS | Γ-centered mesh | k-path |
| LWAVE | .TRUE. | .FALSE. |
| LCHARG | .TRUE. | .FALSE. |

## k-path Explanation

### Γ-M-K-Γ Path (Hexagonal 2D)

```
Γ = (0.0, 0.0, 0.0)    # Brillouin zone center
M = (0.5, 0.0, 0.0)    # Edge midpoint
K = (1/3, 1/3, 0.0)    # Corner
Γ = (0.0, 0.0, 0.0)    # Back to center
```

**Why this path?**
- Standard for hexagonal lattices
- Captures important band features
- VBM usually at K point (TMDC)
- CBM varies (Γ or Q point)

**Number of points:**
- 40 points per segment
- Total: ~120 k-points
- Higher = smoother bands (but slower)

## Expected Results

### Calculation Time:
- **SOC OFF**: 4-6 hours
- Nodes: 2
- Cores: 32

### Output Files:
```
EIGENVAL:  ~500 KB   # Eigenvalues (main result)
OUTCAR:    ~1-2 MB   # Convergence info
DOSCAR:    ~100 MB   # DOS at each k-point
```

### Band Gap:
```
Expected: 1.0-1.5 eV (from DOS)
Type: Likely indirect
VBM: K point (typical TMDC)
CBM: Γ or Q point
```

## Analyzing Results

### Python Script Output:

```bash
python plot_band_structure.py EIGENVAL
```

**Generates:**
- `band_structure.png` - Plot
- Terminal output with:
  - Band gap value
  - Direct/Indirect
  - VBM/CBM k-point locations

### Interpretation:

**Direct Gap:**
```
VBM and CBM at same k-point
→ Efficient optical transitions
→ Good for LEDs, lasers
```

**Indirect Gap:**
```
VBM and CBM at different k-points
→ Requires phonon assistance
→ Lower optical efficiency
→ Typical for TMDCs
```

## Comparison with DOS

### DOS Band Gap:
```
From earlier analysis:
Se-up/S-up: 1.22 eV
```

### Band Structure Gap:
```
Expected: Similar (~1.2 eV)
BUT: More accurate!
```

**Why different?**
- DOS: Approximate (from state density)
- Band: Exact (from E(k))

## Troubleshooting

### Issue: Job fails immediately

**Check:**
```bash
# CHGCAR exists?
ls -lh MoSSe_Se_up__WSSe_S_up/static_scf/CHGCAR

# Should be ~800 MB
# If missing: Re-run static SCF
```

### Issue: No band gap visible

**Possible causes:**
1. Metallic system (unlikely)
2. Fermi level wrong
3. Need wider energy range

**Solution:**
```python
# In plot_band_structure.py, adjust:
ax.set_ylim(-5, 5)  # Instead of (-3, 3)
```

### Issue: Bands look weird

**Check:**
```bash
# Number of k-points
head -6 EIGENVAL | tail -1

# Should be ~120
# If different, check KPOINTS
```

## Next Steps

### After SOC OFF completes:

**Option A: Satisfied with results**
→ Done! Write paper

**Option B: Need SOC**
→ Run SOC ON calculation (next step)

### SOC ON Calculation:

Will provide:
- Spin-orbit splitting
- Valley degeneracy lifting
- More accurate gaps
- Spin textures

**Time: 8-12 hours additional**

## Validation

Before moving to SOC ON:

- [ ] Band structure plot generated?
- [ ] Band gap reasonable (0.5-2 eV)?
- [ ] VBM/CBM locations identified?
- [ ] Consistent with DOS?

## Tips

**For publication:**
- Plot DOS and band structure side-by-side
- Mark VBM/CBM in both
- Compare gap values
- Discuss direct vs indirect

**Common errors:**
- CHGCAR missing → Re-run static SCF with LCHARG=.TRUE.
- Wrong k-path → Check hexagonal symmetry
- Too few k-points → Increase from 40 to 60

## Summary

**What you get:**
- ✅ Accurate band gap
- ✅ Direct/Indirect determination
- ✅ VBM/CBM locations
- ✅ Band dispersion

**What you don't get:**
- ❌ Spin-orbit splitting (need SOC ON)
- ❌ Valley splitting (need SOC ON)
- ❌ Spin textures (need SOC ON)

**Time investment:**
- Setup: 5 minutes
- Calculation: 4-6 hours
- Analysis: 10 minutes

---

**Ready? Run `./setup_band_structure.sh` now!**

After completion, we'll prepare SOC ON calculation.
