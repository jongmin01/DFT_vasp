# PROCAR Analysis for MoSSe/WSSe Janus Heterostructures

Complete analysis toolkit for VASP PROCAR files, specifically designed for TMDC Janus heterostructures.

## Overview

This analysis package provides comprehensive tools to analyze electronic properties from VASP PROCAR files:

### 1. **Orbital-Resolved Band Structure (Fat Bands)**
- Visualize orbital character (s, p, d orbitals) on band structure
- Multi-orbital overlays
- Colormap representations

### 2. **Layer-Resolved Band Structure**
- MoSSe vs WSSe layer contributions
- Band alignment analysis (Type-I vs Type-II)
- Interlayer coupling identification
- VBM/CBM character analysis

### 3. **Combined Band Structure + PDOS**
- Publication-quality combined plots
- Layer-projected DOS
- Comprehensive visualization

### 4. **Valleytronics Analysis**
- K and K' valley identification
- Valley splitting calculations
- Valley degeneracy assessment
- Orbital character at valley points
- Potential for valleytronic applications

### 5. **Band Character Analysis**
- VBM/CBM composition (which layer, which orbital)
- Band gap type (direct/indirect)
- Spatial localization of charge carriers

---

## Quick Start

### Prerequisites

```bash
# Required Python packages
numpy
matplotlib
scipy  # optional, for smoothing
```

### Step 1: Prepare Your Data

Copy VASP output files to analysis directory:

```bash
# Example: Copy files from HPC calculation
cp /path/to/calculation/PROCAR ./90_analysis/MoSSe_Se_up_WSSe_S_up/
cp /path/to/calculation/DOSCAR ./90_analysis/MoSSe_Se_up_WSSe_S_up/
cp /path/to/calculation/CONTCAR ./90_analysis/MoSSe_Se_up_WSSe_S_up/
```

Required files:
- **PROCAR** - Orbital projections
- **DOSCAR** - Density of states (for Fermi energy)
- **POSCAR** or **CONTCAR** - Structure information

### Step 2: Run Comprehensive Analysis

```bash
# Navigate to analysis directory
cd /home/user/DFT_vasp/projects/TMDC_Janus_heterostructure/90_analysis

# Run complete analysis
python scripts/comprehensive_procar_analysis.py \
    --data-dir ./MoSSe_Se_up_WSSe_S_up

# Or specify custom output directory
python scripts/comprehensive_procar_analysis.py \
    --data-dir ./MoSSe_Se_up_WSSe_S_up \
    --output-dir ./analysis_results/0deg_nosoc
```

### Step 3: Review Results

The script generates:

**Plots:**
- `fatband_*.png` - Orbital-resolved band structures
- `layer_resolved_bands.png` - Layer-projected band structure
- `band_dos_combined.png` - Band structure + DOS
- `valley_energies.png` - Valley band energies
- `brillouin_zone_valleys.png` - BZ with K/K' points

**Reports:**
- `COMPREHENSIVE_ANALYSIS_SUMMARY.txt` - Overview of all analyses
- `band_character_analysis.txt` - VBM/CBM character and band alignment
- `valley_analysis_report.txt` - Valleytronics properties

---

## Detailed Usage

### Individual Analysis Modules

You can also run individual analysis modules:

#### 1. Orbital-Resolved Band Structure (Fat Bands)

```bash
# Single orbital
python ../../../../shared/scripts/analysis/plot_fatband.py \
    --procar PROCAR \
    --doscar DOSCAR \
    --orbital dz2 \
    --style single \
    --output fatband_dz2.png

# Multiple orbitals
python ../../../../shared/scripts/analysis/plot_fatband.py \
    --procar PROCAR \
    --doscar DOSCAR \
    --orbitals dz2 dxy dx2 \
    --style multi \
    --output fatband_d_orbitals.png
```

#### 2. Layer-Resolved Analysis

```bash
python ../../../../shared/scripts/analysis/plot_layer_resolved.py \
    --procar PROCAR \
    --poscar CONTCAR \
    --doscar DOSCAR \
    --layers MoSSe WSSe \
    --output layer_bands.png \
    --report band_character.txt
```

#### 3. Band Structure + PDOS

```bash
python ../../../../shared/scripts/analysis/plot_band_pdos.py \
    --procar PROCAR \
    --doscar DOSCAR \
    --output band_pdos.png
```

#### 4. Valleytronics Analysis

```bash
python ../../../../shared/scripts/analysis/analyze_valleytronics.py \
    --procar PROCAR \
    --doscar DOSCAR \
    --output-dir ./valley_analysis
```

---

## Understanding the Results

### Band Character Analysis

The band character analysis reveals:

```
VBM Layer Contributions:
  MoSSe: 0.8234 (82.3%)
  WSSe:  0.1766 (17.7%)

CBM Layer Contributions:
  MoSSe: 0.2145 (21.5%)
  WSSe:  0.7855 (78.5%)
```

**Interpretation:**
- **Type-II band alignment**: VBM on MoSSe layer, CBM on WSSe layer
- Spatial separation of electrons and holes
- Favorable for:
  - Exciton formation
  - Charge separation
  - Photocatalysis applications

### Valley Analysis

For valleytronics applications, key metrics:

```
Valley Splitting:
  VBM: 0.003 eV  ← Nearly degenerate (good!)
  CBM: 0.005 eV  ← Nearly degenerate (good!)
```

**Interpretation:**
- **Small valley splitting** → Good valley degeneracy
- Suitable for valley-based information storage
- Potential for valley Hall effect devices

**For SOC calculations**, valley splitting increases due to spin-valley coupling.

### Orbital Character at Valleys

```
K valley VBM:
  dz2:  0.45 (45%)
  dxy:  0.30 (30%)
  dx2:  0.25 (25%)
```

**Interpretation:**
- VBM primarily d-orbital character
- Dominated by dz2 (out-of-plane orbital)
- Typical for TMDC materials

---

## Analysis Workflow for Different Conditions

### For Twisted Angles (Moiré Patterns)

When analyzing different twist angles:

```bash
# 0 degree
python scripts/comprehensive_procar_analysis.py \
    --data-dir ./0deg_nosoc \
    --output-dir ./analysis/0deg

# 5 degree (when calculated)
python scripts/comprehensive_procar_analysis.py \
    --data-dir ./5deg_nosoc \
    --output-dir ./analysis/5deg
```

Compare:
- Valley splitting vs. twist angle
- Band alignment changes
- Interlayer coupling strength

### For SOC vs. No-SOC

```bash
# Without SOC (current calculation)
python scripts/comprehensive_procar_analysis.py \
    --data-dir ./nosoc \
    --output-dir ./analysis/nosoc

# With SOC (when calculated)
python scripts/comprehensive_procar_analysis.py \
    --data-dir ./soc \
    --output-dir ./analysis/soc
```

Compare:
- Valley splitting (should increase with SOC)
- Spin-valley coupling
- Rashba splitting

---

## Valleytronics Applications

### What to Look For

**1. Valley Degeneracy (for valley memory)**
- VBM/CBM splitting at K vs K' < 10 meV → Excellent
- 10-50 meV → Good
- \> 50 meV → May limit valleytronic performance

**2. Spin-Valley Coupling (requires SOC)**
- Large valley splitting with SOC → Strong spin-valley coupling
- Enables optical valley polarization
- Important for valley-selective excitation

**3. Berry Curvature**
- Opposite at K and K' valleys
- Leads to valley Hall effect
- Can be inferred from orbital character (needs SOC calculation)

### Potential Applications

Based on your results:

1. **Valley-based qubits** - If valley degeneracy is good
2. **Valley photonics** - If optical transitions are valley-selective
3. **Valley Hall effect devices** - If Berry curvature is significant (check with SOC)
4. **Moiré excitons** - If interlayer coupling is tunable with twist angle

---

## Moiré Pattern Analysis

For twisted heterostructures:

### Key Questions to Answer

1. **How does valley splitting change with twist angle?**
   - Plot valley splitting vs. twist angle
   - Identify magic angles

2. **How does band alignment change?**
   - Type-I → Type-II transition?
   - Band offset variation

3. **Interlayer coupling strength**
   - Layer projection at VBM/CBM
   - Hybridization degree

### Analysis Strategy

```bash
# After calculating multiple angles (0°, 5°, 10°, ...)
# Extract valley splitting for each:

for angle in 0 5 10 15 20; do
    python scripts/comprehensive_procar_analysis.py \
        --data-dir ./angle_${angle}deg \
        --output-dir ./analysis/angle_${angle}deg
done

# Then create comparison plots
# (you can write a separate script for this)
```

---

## Troubleshooting

### Common Issues

**1. "PROCAR not found"**
- PROCAR is in .gitignore, so not in git repository
- Copy from HPC calculation results
- Ensure LORBIT = 11 or 12 in INCAR

**2. "Could not identify layers"**
- Check POSCAR/CONTCAR atom ordering
- Manually specify layer atoms if needed
- Edit `identify_tmdc_layers()` function if ordering is unusual

**3. "Valley points not found"**
- Check k-point path includes K and K' points
- For hexagonal: K = (1/3, 1/3, 0), K' = (2/3, 2/3, 0)
- Ensure band structure calculation, not just SCF

**4. "Memory error"**
- PROCAR file very large (especially with many k-points)
- Process on machine with more RAM
- Or modify code to process in chunks

---

## Customization

### Modifying Layer Identification

If automatic layer detection fails, edit `identify_tmdc_layers()` in:
`shared/scripts/analysis/plot_layer_resolved.py`

```python
# Manual layer specification
layers = {
    'MoSSe': [0, 1, 2],      # Mo, S, Se of MoSSe layer
    'WSSe': [3, 4, 5],       # W, S, Se of WSSe layer
}
```

### Adjusting Plot Parameters

In `comprehensive_procar_analysis.py`, you can adjust:
- Energy range: `energy_range=(-3, 3)` → `(-5, 5)` for wider view
- Colormap: `cmap='hot'` → `'viridis'`, `'plasma'`, etc.
- DPI: `dpi=300` → `600` for higher resolution

---

## Output Files Summary

### Automatically Generated

| File | Description |
|------|-------------|
| `COMPREHENSIVE_ANALYSIS_SUMMARY.txt` | Overview of all analyses |
| `fatband_dz2.png` | dz2 orbital fat band |
| `fatband_dxy.png` | dxy orbital fat band |
| `fatband_d_orbitals.png` | Multi-orbital d-orbitals |
| `fatband_d_total_colormap.png` | Total d-character colormap |
| `layer_resolved_bands.png` | Layer-projected band structure |
| `band_character_analysis.txt` | VBM/CBM character, Type-I/II |
| `band_dos_combined.png` | Band + DOS combined plot |
| `valley_energies.png` | Valley band energies |
| `brillouin_zone_valleys.png` | BZ with marked valleys |
| `valley_analysis_report.txt` | Detailed valley properties |

---

## Next Steps

### Recommended Analysis Sequence

1. **Initial characterization** (0° structure, no SOC) ✓ You are here
   - Run comprehensive analysis
   - Identify band alignment type
   - Check valley degeneracy

2. **SOC calculation** (if valleytronics is important)
   - Run VASP with LSORBIT = .TRUE.
   - Re-run analysis
   - Compare valley splitting

3. **Twisted angle series** (for Moiré physics)
   - Calculate 5°, 10°, 15°, 20°, ...
   - Analyze each with this toolkit
   - Plot trends vs. twist angle

4. **Publication plots**
   - Use generated plots as starting point
   - Customize for publication quality
   - Combine multiple analyses

---

## Support and Development

### Module Structure

```
shared/scripts/analysis/
├── procar_parser.py           # PROCAR file parser
├── plot_fatband.py            # Orbital-resolved bands
├── plot_layer_resolved.py     # Layer-resolved bands
├── plot_band_pdos.py          # Band + PDOS combined
└── analyze_valleytronics.py   # Valley physics analysis

projects/TMDC_Janus_heterostructure/90_analysis/
├── scripts/
│   └── comprehensive_procar_analysis.py  # Main analysis script
└── README_PROCAR_ANALYSIS.md             # This file
```

### Adding New Features

To add new analysis:
1. Create new module in `shared/scripts/analysis/`
2. Import in `comprehensive_procar_analysis.py`
3. Add new method to `ComprehensiveProcarAnalysis` class
4. Call in `run_all_analyses()`

---

## References

For background on TMDC valleytronics and Janus structures:

1. **Valley physics in TMDCs:**
   - Xiao et al., PRL 108, 196802 (2012)
   - Mak et al., Nature Nanotech. 7, 494 (2012)

2. **Janus TMDC structures:**
   - Lu et al., Nature Nanotech. 12, 744 (2017)
   - Zhang et al., ACS Nano 11, 8192 (2017)

3. **Moiré patterns in heterostructures:**
   - Wu et al., Nature 567, 323 (2019)
   - Tran et al., Nature 567, 71 (2019)

---

## Version History

- **v1.0** (2025-01-23): Initial release
  - Complete PROCAR analysis toolkit
  - Orbital, layer, and valley analyses
  - Comprehensive reporting

---

Good luck with your analysis! 🚀

For questions or issues, check the generated report files first - they contain detailed interpretations.
