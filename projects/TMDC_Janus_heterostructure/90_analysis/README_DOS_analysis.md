# DOS Analysis Guide

## Files Generated

1. **[plot_dos.py](computer:///mnt/user-data/outputs/plot_dos.py)** - Single structure DOS plotting
2. **[compare_all_dos.py](computer:///mnt/user-data/outputs/compare_all_dos.py)** - Compare all 4 structures

## Quick Start

### Method 1: Analyze Single Structure

```bash
# Move to calculation directory
cd /path/to/your/calculations

# Plot DOS for one structure
python plot_dos.py MoSSe_Se_up__WSSe_S_up/static_scf/DOSCAR
```

**Output:**
- `DOS.png` - DOS plot
- Band gap estimate printed to terminal

### Method 2: Compare All Structures (Recommended)

```bash
# In parent directory with all MoSSe_* folders
python compare_all_dos.py
```

**Output:**
- `DOS_all_structures.png` - 2×2 grid of individual plots
- `DOS_comparison.png` - Overlay comparison
- `DOS_analysis_summary.txt` - Numerical results

## What You'll Get

### 1. DOS Plots

Each plot shows:
- **Blue curve**: Density of states
- **Black dashed line**: Fermi level (E_F = 0)
- **Orange dotted lines**: VBM and CBM (if gap exists)
- **Gray shaded area**: Band gap region
- **Text box**: Estimated band gap value

### 2. Band Gap Estimation

The script estimates:
- **Band gap**: Energy difference between VBM and CBM
- **VBM**: Valence band maximum
- **CBM**: Conduction band minimum

**Note:** This is an **estimate** from DOS. For accurate determination, band structure is needed.

### 3. Comparison Plots

**Individual plots (2×2):**
- Side-by-side comparison of all stackings
- Easy to see differences

**Overlay plot:**
- All DOS on same axes
- Direct comparison of features

## Understanding the Results

### Band Gap Interpretation

```
Gap < 0.5 eV  : Narrow gap semiconductor
Gap 0.5-2 eV  : Moderate gap (typical for TMDCs)
Gap 2-5 eV    : Wide gap semiconductor
Gap > 5 eV    : Very wide gap (check if correct)
```

### Fermi Level Position

```
E_F near VBM  : p-type character
E_F at midgap : Intrinsic semiconductor
E_F near CBM  : n-type character
```

### DOS Shape

**Near E_F:**
- Sharp peaks = Localized states
- Broad features = Delocalized states

**d-orbital character (expected for Mo/W):**
- Multiple peaks between -2 to 2 eV

## Expected Results for TMDC Heterostructures

### Typical values:
```
Band gap: 1.0-2.0 eV (indirect)
DOS peaks: Multiple peaks from Mo/W d-orbitals
Fermi level: Slightly negative (n-type tendency)
```

### Your structures:
Based on relaxation energies, expect:
- All 4 structures: Similar band gaps (within 0.1-0.2 eV)
- Structure 2 (most stable): Possibly slightly different DOS shape
- Overall: Semiconductor behavior

## Advanced Analysis (Optional)

### 1. Extract Projected DOS (PDOS)

If you need orbital-resolved DOS:

```bash
# DOSCAR has PDOS if LORBIT=11 was set
# Parse columns 2-end for different orbitals
```

### 2. Integration

Calculate number of states:

```python
from scipy.integrate import trapz

# Integrate DOS up to Fermi level
n_electrons = trapz(dos[energy < 0], energy[energy < 0])
```

### 3. Effective Mass Estimation

Requires band structure, not DOS.

## Troubleshooting

### Issue: No plot generated

**Check:**
```bash
# 1. Python and matplotlib installed?
python --version
python -c "import matplotlib"

# 2. DOSCAR exists?
ls -lh MoSSe_*/static_scf/DOSCAR

# 3. Run from correct directory?
pwd  # Should be parent of MoSSe_* directories
```

### Issue: Band gap = 0 or negative

**Possible causes:**
1. System is metallic (check relaxation results)
2. DOS threshold too high (adjust in script)
3. Fermi level determination issue

**Solution:**
```python
# In script, adjust threshold:
gap, vbm, cbm = estimate_band_gap(energy, dos, threshold=0.005)
```

### Issue: Unrealistic band gap (> 5 eV)

**Check:**
1. DOSCAR file correct format?
2. Energy range reasonable?
3. Static SCF converged properly?

## Python Requirements

### Minimal:
```bash
pip install numpy matplotlib
```

### Full (recommended):
```bash
pip install numpy matplotlib scipy
```

## Alternative Tools

If Python not available:

### 1. p4vasp
```bash
p4v MoSSe_Se_up__WSSe_S_up/static_scf/vasprun.xml
# GUI: Electronic -> DOS
```

### 2. vaspkit
```bash
cd MoSSe_Se_up__WSSe_S_up/static_scf
vaspkit
# Select: 11 (Density of States)
```

### 3. VESTA
```bash
# Can visualize DOS from vasprun.xml
```

## Output Files Summary

After running scripts:

```
your_directory/
├── MoSSe_Se_up__WSSe_Se_up/
│   └── static_scf/
│       ├── DOSCAR (input)
│       └── DOS.png (output, if individual)
├── MoSSe_Se_up__WSSe_S_up/
│   └── static_scf/DOSCAR
├── MoSSe_S_up__WSSe_Se_up/
│   └── static_scf/DOSCAR
├── MoSSe_S_up__WSSe_S_up/
│   └── static_scf/DOSCAR
├── DOS_all_structures.png (2×2 grid)
├── DOS_comparison.png (overlay)
└── DOS_analysis_summary.txt (numerical)
```

## Next Steps

After DOS analysis:

1. **Interpret results** - Band gaps, Fermi level
2. **Compare stackings** - Which is most stable?
3. **Decide**: Need band structure?
   - Yes: For E(k), effective mass, valley
   - No: DOS sufficient for basic analysis
4. **Decide**: Need SOC?
   - Yes: For spin-orbit splitting
   - No: SOC OFF sufficient

## Example Workflow

```bash
# 1. Upload scripts to server
cd /home/jbaek/vasp/hetero/Janus/relax_0deg_structures

# 2. Run comparison
python compare_all_dos.py

# 3. Download plots
# Use scp or file transfer to get PNG files

# 4. Analyze
# Open images, read summary.txt

# 5. Decide next steps
# Band structure? SOC? Or finish here?
```

## Tips

- **Resolution**: Plots are 300 DPI (publication quality)
- **Customization**: Edit Python scripts for different ranges, colors
- **Energy range**: Default -5 to +5 eV (change in script if needed)
- **DOS threshold**: Adjust if gap detection fails

---

**Ready to start? Upload scripts and run `compare_all_dos.py`!**
