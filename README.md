# VASP DFT Calculations Repository

This repository contains DFT (Density Functional Theory) calculations performed using VASP, organized by research projects and materials.

## Repository Structure

```
DFT_vasp/
├── projects/              # Research projects organized by topic
│   ├── MoS2_monolayer/   # MoS2 monolayer calculations
│   └── TMDC_Janus_heterostructure/  # Janus TMDC heterostructure studies
│       └── twisted_angle/           # Twisted angle research
├── scripts/              # Analysis and utility scripts
│   ├── structure_generation/  # Structure creation scripts
│   ├── analysis/             # Data analysis scripts
│   └── plotting/             # Visualization scripts
└── docs/                 # Documentation and notes
```

## Workflow

1. **Structure Generation**: Create initial structures locally using scripts in `scripts/structure_generation/`
2. **HPC Calculation**: Transfer input files to ANL carbon HPC for VASP calculations
3. **Data Transfer**: Download essential files (INCAR, POSCAR, KPOINTS, CONTCAR, vasprun.xml) to local repository
4. **Analysis & Plotting**: Analyze results and create plots locally using scripts

## Calculation Types

Each material/system typically includes the following calculation directories:

- `relax/` - Structural optimization (relaxation)
- `scf/` - Self-consistent field calculation
- `band/` - Band structure calculation
- `dos/` - Density of states calculation

## Essential Files Kept in Repository

- **INCAR** - Input parameters for VASP
- **POSCAR** - Initial atomic structure
- **KPOINTS** - K-point mesh settings
- **CONTCAR** - Relaxed structure (contains material information)
- **vasprun.xml** - Main output file with calculation results
- **EIGENVAL** - Eigenvalues for band structure

## Files Excluded (in .gitignore)

Large files and licensed files are excluded from the repository:
- **POTCAR** - Pseudopotentials (excluded due to VASP license restrictions for public repositories)
- **OUTCAR** - Large output file (data is in vasprun.xml)
- **WAVECAR** - Very large wavefunction files
- **CHG/CHGCAR** - Large charge density files
- **DOSCAR** - DOS data (data is in vasprun.xml)

## Current Projects

### 1. MoS2 Monolayer
**Status**: Completed
- Relaxation ✓
- SCF ✓
- Band structure ✓
- DOS ✓

### 2. TMDC Janus Heterostructure - Twisted Angle Study
**Status**: In progress
- Research focus: Effect of twist angles on electronic properties and Moiré pattern formation
- Systems: MoSSe/WSSe heterostructures
- Initial test angles: 0°, 5°, 10°, 15°, 20°
- Current: 0° relaxation in progress
- Future plans:
  - Extend to additional angles based on Moiré pattern analysis
  - Investigate other TMDC Janus combinations
  - Study other semiconductor heterostructures

## Adding New Projects

1. Create project directory under `projects/`
2. Organize by research topic/material system
3. Create subdirectories for different calculations (relax, scf, band, dos)
4. Add project description to this README

## Notes

- Calculations performed on ANL carbon HPC cluster
- VASP version: [Add your VASP version]
- For questions or collaboration, contact [Your contact info]
