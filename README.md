# VASP DFT Calculations Repository

This repository contains DFT (Density Functional Theory) calculations performed using VASP, organized by research projects and materials.

## Repository Structure

```
vasp/
├── projects/                       # One directory per research project
│   ├── MoS2_monolayer/
│   │   ├── 00_structure_generation/  # Seeds + project-specific generators
│   │   ├── 10_relax/                # All relaxation runs
│   │   ├── 20_static_scf/           # Static SCF calculations
│   │   ├── 30_band/                 # Band-structure workflows
│   │   ├── 40_dos/                  # DOS workflows
│   │   └── 90_analysis/             # Plots, notebooks, summaries
│   └── TMDC_Janus_heterostructure/  # Shares the same stage layout
├── shared/
│   ├── scripts/                     # Reusable helpers (HPC, plotting, builders)
│   ├── templates/                   # Canonical INCAR/KPOINTS/job files
│   └── tools/                       # Cross-project utilities (visualization, etc.)
├── structures/
│   ├── seed_structures/             # Global library of initial POSCARs
│   └── reference/JARVIS/            # Reference data pulled from JARVIS
└── docs/                            # Documentation and workflow notes
```

## Workflow

1. **Structure Generation**: Create initial structures locally using scripts in `projects/<project>/00_structure_generation/` or reusable helpers in `shared/scripts/structure_generation/`
2. **HPC Calculation**: Transfer input files to ANL carbon HPC for VASP calculations
3. **Data Transfer**: Download essential files (INCAR, POSCAR, KPOINTS, CONTCAR, vasprun.xml) to local repository
4. **Analysis & Plotting**: Analyze results and create plots locally using scripts

## Calculation Types

Each project now follows a consistent numbered pipeline so you can immediately locate files:

1. `00_structure_generation/` - Structural seeds, builders, and setup scripts
2. `10_relax/` - Structural optimization (relaxations)
3. `20_static_scf/` - Self-consistent field calculations (static runs)
4. `30_band/` - Band structure calculations and post-processing
5. `40_dos/` - Density of states calculations
6. `90_analysis/` - Project-specific analysis, plots, and summaries

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
- **vasprun.xml** (selected cases) - Extremely large `vasprun.xml` snapshots under  
  `projects/MoS2_monolayer/90_analysis/MoS2_strain_plots/*/soc/` and  
  `projects/TMDC_Janus_heterostructure/20_static_scf/*/scf/` exceed GitHub's 100 MB limit,  
  so they are stored locally only (see `.gitignore` for exact paths).  Derived plots and
  summaries remain available in `90_analysis/`.

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
3. Create the standard stage folders (`00_structure_generation`, `10_relax`, `20_static_scf`, `30_band`, `40_dos`, `90_analysis`)
4. Add project description to this README

## Notes

- Calculations performed on ANL carbon HPC cluster
- VASP version: [Add your VASP version]
- For questions or collaboration, contact [Your contact info]
