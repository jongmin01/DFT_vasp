#!/usr/bin/env python3
"""
comprehensive_procar_analysis.py - Complete PROCAR analysis for MoSSe/WSSe heterostructures

This script performs ALL available PROCAR analyses:
1. Orbital-resolved band structure (fat bands)
2. Layer-resolved band structure
3. Combined band structure + PDOS
4. Band character analysis (VBM/CBM)
5. Valleytronics analysis (K/K' valleys)
6. Comprehensive report generation

Usage:
    python comprehensive_procar_analysis.py --data-dir <path-to-vasp-outputs>

Expected files in data-dir:
    - PROCAR
    - DOSCAR
    - POSCAR (or CONTCAR)
    - EIGENVAL (optional)
    - OUTCAR (optional)
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
import os
from pathlib import Path
import argparse
from typing import Dict, List, Optional

# Add shared scripts to path
script_dir = Path(__file__).parent.parent.parent.parent / 'shared' / 'scripts' / 'analysis'
sys.path.insert(0, str(script_dir))

# Import analysis modules
import procar_parser as pp
from plot_fatband import (
    plot_fatband_single_orbital,
    plot_fatband_multiorbital,
    plot_fatband_colormap,
    read_fermi_energy
)
from plot_layer_resolved import (
    identify_tmdc_layers,
    plot_layer_resolved_bands,
    analyze_band_character
)
from plot_band_pdos import (
    read_doscar_total,
    plot_band_and_pdos,
    plot_band_with_layer_pdos
)
from analyze_valleytronics import (
    identify_valley_points,
    analyze_valley_band_structure,
    compute_valley_splitting,
    plot_valley_band_energies,
    plot_brillouin_zone_with_valleys,
    analyze_orbital_character_at_valleys,
    generate_valley_analysis_report
)


class ComprehensiveProcarAnalysis:
    """Complete PROCAR analysis pipeline"""

    def __init__(self, data_dir: str, output_dir: str = None):
        """
        Initialize analysis

        Args:
            data_dir: Directory containing VASP output files
            output_dir: Directory for output files (default: data_dir/procar_analysis)
        """
        self.data_dir = Path(data_dir)
        self.output_dir = Path(output_dir) if output_dir else self.data_dir / 'procar_analysis'

        # Create output directory
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # File paths
        self.procar_file = self.data_dir / 'PROCAR'
        self.doscar_file = self.data_dir / 'DOSCAR'
        self.poscar_file = self._find_structure_file()

        # Check required files
        self._check_required_files()

        # Data containers
        self.procar_data = None
        self.efermi = None
        self.layer_atoms = None
        self.valley_indices = None

        print(f"Comprehensive PROCAR Analysis")
        print(f"=" * 80)
        print(f"Data directory: {self.data_dir}")
        print(f"Output directory: {self.output_dir}")
        print(f"=" * 80)

    def _find_structure_file(self) -> Path:
        """Find POSCAR or CONTCAR file"""
        if (self.data_dir / 'CONTCAR').exists():
            return self.data_dir / 'CONTCAR'
        elif (self.data_dir / 'POSCAR').exists():
            return self.data_dir / 'POSCAR'
        else:
            raise FileNotFoundError("Neither POSCAR nor CONTCAR found")

    def _check_required_files(self):
        """Check if required files exist"""
        required = {
            'PROCAR': self.procar_file,
            'DOSCAR': self.doscar_file,
            'POSCAR/CONTCAR': self.poscar_file
        }

        missing = []
        for name, path in required.items():
            if not path.exists():
                missing.append(name)

        if missing:
            raise FileNotFoundError(f"Missing required files: {', '.join(missing)}")

        print("✓ All required files found")

    def parse_data(self):
        """Parse all input data"""
        print(f"\n{'='*80}")
        print("STEP 1: Parsing input files")
        print(f"{'='*80}")

        # Parse PROCAR
        print(f"\n1.1 Parsing PROCAR: {self.procar_file}")
        parser = pp.ProcarParser(str(self.procar_file))
        self.procar_data = parser.parse()

        print(f"  ✓ K-points: {self.procar_data.nkpts}")
        print(f"  ✓ Bands: {self.procar_data.nbands}")
        print(f"  ✓ Ions: {self.procar_data.nions}")
        print(f"  ✓ Orbitals: {self.procar_data.orbital_names}")
        print(f"  ✓ Spin-polarized: {self.procar_data.is_spin_polarized}")
        print(f"  ✓ SOC: {self.procar_data.is_soc}")

        # Read Fermi energy
        print(f"\n1.2 Reading Fermi energy from: {self.doscar_file}")
        self.efermi = read_fermi_energy(str(self.doscar_file))
        print(f"  ✓ Fermi energy: {self.efermi:.4f} eV")

        # Identify layers
        print(f"\n1.3 Identifying layers from: {self.poscar_file}")
        self.layer_atoms = identify_tmdc_layers(str(self.poscar_file))

        # Identify valley points
        print(f"\n1.4 Identifying valley points (K, K')")
        self.valley_indices = identify_valley_points(self.procar_data.kpoints)

        for valley_name, idx in self.valley_indices.items():
            if not valley_name.endswith('_coord'):
                print(f"  ✓ {valley_name}: k-point index {idx}")

    def run_orbital_analysis(self):
        """Run orbital-resolved band structure analysis"""
        print(f"\n{'='*80}")
        print("STEP 2: Orbital-resolved band structure (Fat bands)")
        print(f"{'='*80}")

        energies = (self.procar_data.energies[:, :, 0] if self.procar_data.is_spin_polarized
                   else self.procar_data.energies)

        # Define orbital groups
        d_orbitals = ['dxy', 'dyz', 'dz2', 'dxz', 'dx2']
        p_orbitals = ['py', 'pz', 'px']

        # Available orbitals
        available_d = [orb for orb in d_orbitals if orb in self.procar_data.orbital_names]
        available_p = [orb for orb in p_orbitals if orb in self.procar_data.orbital_names]

        # 2.1: Individual d-orbitals (most important for TMDC)
        if available_d:
            print(f"\n2.1 Plotting d-orbital fat bands...")

            for orb in available_d[:3]:  # Plot first 3 d-orbitals
                proj = self._get_orbital_projection_summed(orb)

                plot_fatband_single_orbital(
                    self.procar_data.kpoints,
                    energies,
                    proj,
                    orb,
                    self.efermi,
                    title=f'Fat Band Structure - {orb} orbital',
                    output=str(self.output_dir / f'fatband_{orb}.png')
                )
                print(f"  ✓ {orb} fat band plot saved")

        # 2.2: Multi-orbital plot (all d-orbitals)
        if available_d:
            print(f"\n2.2 Plotting multi-orbital fat bands (d-orbitals)...")

            proj_dict = {}
            for orb in available_d:
                proj_dict[orb] = self._get_orbital_projection_summed(orb)

            plot_fatband_multiorbital(
                self.procar_data.kpoints,
                energies,
                proj_dict,
                self.efermi,
                title='Fat Band Structure - d-orbitals',
                output=str(self.output_dir / 'fatband_d_orbitals.png')
            )
            print(f"  ✓ Multi-orbital d-orbital fat band plot saved")

        # 2.3: Colormap plot for total d-orbital character
        if available_d:
            print(f"\n2.3 Plotting colormap fat bands (total d-character)...")

            # Sum all d-orbital projections
            d_proj_total = np.zeros_like(energies)
            for orb in available_d:
                d_proj_total += self._get_orbital_projection_summed(orb)

            plot_fatband_colormap(
                self.procar_data.kpoints,
                energies,
                d_proj_total,
                "d-orbital",
                self.efermi,
                title='Fat Band Structure - Total d-orbital Character',
                output=str(self.output_dir / 'fatband_d_total_colormap.png')
            )
            print(f"  ✓ d-orbital colormap plot saved")

    def run_layer_analysis(self):
        """Run layer-resolved band structure analysis"""
        print(f"\n{'='*80}")
        print("STEP 3: Layer-resolved band structure")
        print(f"{'='*80}")

        energies = (self.procar_data.energies[:, :, 0] if self.procar_data.is_spin_polarized
                   else self.procar_data.energies)

        # Get layer projections
        layer_projections = {}
        for layer_name in ['MoSSe', 'WSSe']:
            if layer_name in self.layer_atoms:
                print(f"\n3.1 Computing {layer_name} layer projection...")
                proj = self._get_atom_projection(self.layer_atoms[layer_name])
                layer_projections[layer_name] = proj
                print(f"  ✓ {layer_name} projection computed")

        # Plot layer-resolved bands
        if layer_projections:
            print(f"\n3.2 Plotting layer-resolved band structure...")

            plot_layer_resolved_bands(
                self.procar_data.kpoints,
                energies,
                layer_projections,
                self.efermi,
                title='Layer-Resolved Band Structure - MoSSe/WSSe',
                output=str(self.output_dir / 'layer_resolved_bands.png')
            )
            print(f"  ✓ Layer-resolved band structure saved")

            # Band character analysis
            print(f"\n3.3 Analyzing band character (VBM/CBM)...")

            analyze_band_character(
                energies,
                layer_projections,
                self.efermi,
                output_file=str(self.output_dir / 'band_character_analysis.txt')
            )
            print(f"  ✓ Band character analysis saved")

    def run_band_pdos_analysis(self):
        """Run combined band structure + PDOS analysis"""
        print(f"\n{'='*80}")
        print("STEP 4: Combined band structure + PDOS")
        print(f"{'='*80}")

        energies = (self.procar_data.energies[:, :, 0] if self.procar_data.is_spin_polarized
                   else self.procar_data.energies)

        # Read DOSCAR
        print(f"\n4.1 Reading total DOS from DOSCAR...")
        dos_energy, dos_total, _ = read_doscar_total(str(self.doscar_file))
        print(f"  ✓ DOS data loaded")

        # Plot band + total DOS
        print(f"\n4.2 Plotting band structure + total DOS...")

        plot_band_and_pdos(
            self.procar_data.kpoints,
            energies,
            dos_energy,
            {'Total DOS': dos_total},
            self.efermi,
            title='Band Structure and DOS - MoSSe/WSSe',
            output=str(self.output_dir / 'band_dos_combined.png')
        )
        print(f"  ✓ Combined band+DOS plot saved")

        # Plot band + layer-PDOS (if we can compute PDOS from PROCAR)
        # Note: This would require more sophisticated PDOS calculation
        # For now, we skip this part

    def run_valleytronics_analysis(self):
        """Run valleytronics analysis"""
        print(f"\n{'='*80}")
        print("STEP 5: Valleytronics analysis (K/K' valleys)")
        print(f"{'='*80}")

        energies = (self.procar_data.energies[:, :, 0] if self.procar_data.is_spin_polarized
                   else self.procar_data.energies)

        # Analyze valley band structure
        print(f"\n5.1 Analyzing valley band structure...")
        valley_results = analyze_valley_band_structure(
            energies,
            None,  # No orbital projections for now
            self.valley_indices,
            self.efermi
        )
        print(f"  ✓ Valley band structure analyzed")

        # Compute valley splitting
        print(f"\n5.2 Computing valley splitting...")
        splitting = compute_valley_splitting(valley_results)

        if splitting:
            if 'vbm_splitting' in splitting:
                print(f"  ✓ VBM valley splitting: {splitting['vbm_splitting']:.4f} eV")
            if 'cbm_splitting' in splitting:
                print(f"  ✓ CBM valley splitting: {splitting['cbm_splitting']:.4f} eV")

        # Analyze orbital character at valleys
        print(f"\n5.3 Analyzing orbital character at K/K' valleys...")
        orbital_analysis = analyze_orbital_character_at_valleys(
            self.procar_data,
            self.valley_indices
        )
        print(f"  ✓ Orbital character at valleys analyzed")

        # Plot valley energies
        print(f"\n5.4 Plotting valley band energies...")
        plot_valley_band_energies(
            valley_results,
            output=str(self.output_dir / 'valley_energies.png'),
            title='Valley Band Energies - MoSSe/WSSe'
        )
        print(f"  ✓ Valley energy plot saved")

        # Plot Brillouin zone
        print(f"\n5.5 Plotting Brillouin zone with valleys...")
        plot_brillouin_zone_with_valleys(
            self.valley_indices,
            output=str(self.output_dir / 'brillouin_zone_valleys.png')
        )
        print(f"  ✓ Brillouin zone plot saved")

        # Generate valley analysis report
        print(f"\n5.6 Generating valley analysis report...")
        generate_valley_analysis_report(
            valley_results,
            splitting,
            orbital_analysis,
            output_file=str(self.output_dir / 'valley_analysis_report.txt')
        )

    def generate_summary_report(self):
        """Generate comprehensive summary report"""
        print(f"\n{'='*80}")
        print("STEP 6: Generating comprehensive summary report")
        print(f"{'='*80}")

        report = []

        report.append("=" * 80)
        report.append("COMPREHENSIVE PROCAR ANALYSIS REPORT")
        report.append("MoSSe/WSSe Janus TMDC Heterostructure")
        report.append("=" * 80)

        # System information
        report.append("\n" + "=" * 80)
        report.append("1. SYSTEM INFORMATION")
        report.append("=" * 80)
        report.append(f"\nData directory: {self.data_dir}")
        report.append(f"Output directory: {self.output_dir}")
        report.append(f"\nNumber of k-points: {self.procar_data.nkpts}")
        report.append(f"Number of bands: {self.procar_data.nbands}")
        report.append(f"Number of ions: {self.procar_data.nions}")
        report.append(f"Orbitals included: {', '.join(self.procar_data.orbital_names)}")
        report.append(f"Spin-polarized: {self.procar_data.is_spin_polarized}")
        report.append(f"SOC included: {self.procar_data.is_soc}")
        report.append(f"Fermi energy: {self.efermi:.4f} eV")

        # Layer information
        report.append("\n" + "=" * 80)
        report.append("2. LAYER COMPOSITION")
        report.append("=" * 80)

        for layer_name, atom_indices in self.layer_atoms.items():
            if len(atom_indices) <= 10:  # Only show main layers
                report.append(f"\n{layer_name}:")
                report.append(f"  Atom indices: {atom_indices}")
                report.append(f"  Number of atoms: {len(atom_indices)}")

        # Analysis outputs
        report.append("\n" + "=" * 80)
        report.append("3. GENERATED OUTPUT FILES")
        report.append("=" * 80)

        report.append("\nOrbital-resolved analysis:")
        report.append("  - fatband_*.png: Individual orbital fat band plots")
        report.append("  - fatband_d_orbitals.png: Multi-orbital d-orbital plot")
        report.append("  - fatband_d_total_colormap.png: Total d-character colormap")

        report.append("\nLayer-resolved analysis:")
        report.append("  - layer_resolved_bands.png: Layer-projected band structure")
        report.append("  - band_character_analysis.txt: VBM/CBM character analysis")

        report.append("\nCombined plots:")
        report.append("  - band_dos_combined.png: Band structure + DOS")

        report.append("\nValleytronics analysis:")
        report.append("  - valley_energies.png: Valley band energies")
        report.append("  - brillouin_zone_valleys.png: Brillouin zone with K/K' points")
        report.append("  - valley_analysis_report.txt: Detailed valley analysis")

        report.append("\n" + "=" * 80)
        report.append("4. KEY FINDINGS")
        report.append("=" * 80)

        report.append("\nFor detailed results, see individual analysis files:")
        report.append("  - band_character_analysis.txt: Band alignment (Type-I/Type-II)")
        report.append("  - valley_analysis_report.txt: Valleytronics potential")

        report.append("\n" + "=" * 80)

        # Write report
        report_file = self.output_dir / 'COMPREHENSIVE_ANALYSIS_SUMMARY.txt'
        with open(report_file, 'w') as f:
            f.write('\n'.join(report))

        print(f"\n✓ Summary report saved to: {report_file}")
        print("\n" + '\n'.join(report))

    def _get_orbital_projection_summed(self, orbital_name: str) -> np.ndarray:
        """Get orbital projection summed over atoms"""
        parser = pp.ProcarParser(str(self.procar_file))
        parser.data = self.procar_data

        proj = parser.get_orbital_projection(orbital_name)

        if self.procar_data.is_spin_polarized:
            proj = proj[:, :, :, 0]  # Use spin-up

        # Sum over atoms
        proj_sum = proj.sum(axis=2)

        return proj_sum

    def _get_atom_projection(self, atom_indices: List[int]) -> np.ndarray:
        """Get atom projection summed over orbitals"""
        parser = pp.ProcarParser(str(self.procar_file))
        parser.data = self.procar_data

        proj = parser.get_atom_projection(atom_indices)

        if self.procar_data.is_spin_polarized:
            proj = proj[:, :, 0]  # Use spin-up

        return proj

    def run_all_analyses(self):
        """Run complete analysis pipeline"""
        print("\n" + "=" * 80)
        print("STARTING COMPREHENSIVE PROCAR ANALYSIS")
        print("=" * 80)

        # Step 1: Parse data
        self.parse_data()

        # Step 2: Orbital analysis
        self.run_orbital_analysis()

        # Step 3: Layer analysis
        self.run_layer_analysis()

        # Step 4: Band + PDOS
        self.run_band_pdos_analysis()

        # Step 5: Valleytronics
        self.run_valleytronics_analysis()

        # Step 6: Summary report
        self.generate_summary_report()

        print("\n" + "=" * 80)
        print("COMPREHENSIVE ANALYSIS COMPLETE!")
        print("=" * 80)
        print(f"\nAll outputs saved to: {self.output_dir}")
        print("\nNext steps:")
        print("  1. Review COMPREHENSIVE_ANALYSIS_SUMMARY.txt for overview")
        print("  2. Check individual analysis reports for detailed results")
        print("  3. Examine generated plots for visualization")


def main():
    """Main function"""
    parser = argparse.ArgumentParser(
        description='Comprehensive PROCAR analysis for MoSSe/WSSe heterostructures',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Analyze VASP outputs in current directory
  python comprehensive_procar_analysis.py --data-dir .

  # Analyze outputs in specific directory
  python comprehensive_procar_analysis.py --data-dir /path/to/vasp/outputs

  # Specify custom output directory
  python comprehensive_procar_analysis.py --data-dir . --output-dir ./analysis_results
        """
    )

    parser.add_argument('--data-dir', required=True,
                       help='Directory containing VASP output files (PROCAR, DOSCAR, POSCAR)')
    parser.add_argument('--output-dir', default=None,
                       help='Output directory (default: data-dir/procar_analysis)')

    args = parser.parse_args()

    # Create analysis object
    analysis = ComprehensiveProcarAnalysis(
        data_dir=args.data_dir,
        output_dir=args.output_dir
    )

    # Run all analyses
    analysis.run_all_analyses()


if __name__ == '__main__':
    main()
