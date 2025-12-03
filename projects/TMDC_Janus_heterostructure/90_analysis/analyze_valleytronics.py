#!/usr/bin/env python3
"""
analyze_valleytronics.py - Valley physics analysis for TMDC heterostructures

Analyzes valley-related properties in TMDC materials:
- K and K' valley identification and properties
- Valley splitting and degeneracy
- Orbital character at K/K' points (d orbitals)
- Berry curvature indicators
- Valley polarization analysis
- Spin-valley coupling (if SOC is included)

For MoSSe/WSSe Janus heterostructures
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch
import sys
import os
from typing import List, Dict, Optional, Tuple

# Import local modules
import procar_parser as pp


def identify_valley_points(kpoints: np.ndarray,
                          lattice_constant: float = 3.18) -> Dict[str, int]:
    """
    Identify K and K' valley points in the Brillouin zone

    For hexagonal lattice:
    K  = (1/3, 1/3, 0) in fractional coordinates
    K' = (2/3, 2/3, 0) = (-1/3, -1/3, 0) in fractional coordinates

    Args:
        kpoints: k-point coordinates [nkpts, 3]
        lattice_constant: Lattice constant (angstrom)

    Returns:
        Dictionary with K and K' indices and Gamma point
    """
    # Define high-symmetry points in fractional coordinates
    K_frac = np.array([1/3, 1/3, 0])
    Kprime_frac = np.array([2/3, 2/3, 0])  # or [-1/3, -1/3, 0]
    Gamma_frac = np.array([0, 0, 0])

    # Also check equivalent points
    K_equiv = [
        np.array([1/3, 1/3, 0]),
        np.array([2/3, -1/3, 0]),
        np.array([-1/3, 2/3, 0]),
    ]

    Kprime_equiv = [
        np.array([2/3, 2/3, 0]),
        np.array([-1/3, -2/3, 0]),
        np.array([-2/3, 1/3, 0]),
    ]

    valley_points = {}
    tolerance = 0.05  # Tolerance for k-point matching

    # Find Gamma point
    gamma_dist = np.linalg.norm(kpoints - Gamma_frac, axis=1)
    gamma_idx = np.argmin(gamma_dist)
    if gamma_dist[gamma_idx] < tolerance:
        valley_points['Gamma'] = gamma_idx
        valley_points['Gamma_coord'] = kpoints[gamma_idx]

    # Find K point
    k_distances = []
    for k_equiv in K_equiv:
        distances = np.linalg.norm(kpoints - k_equiv, axis=1)
        k_distances.extend(distances)

    k_idx = np.argmin(k_distances)
    k_idx_actual = k_idx % len(kpoints)

    min_dist = np.min(k_distances)
    if min_dist < tolerance:
        valley_points['K'] = k_idx_actual
        valley_points['K_coord'] = kpoints[k_idx_actual]

    # Find K' point
    kprime_distances = []
    for kprime_equiv in Kprime_equiv:
        distances = np.linalg.norm(kpoints - kprime_equiv, axis=1)
        kprime_distances.extend(distances)

    kprime_idx = np.argmin(kprime_distances)
    kprime_idx_actual = kprime_idx % len(kpoints)

    min_dist_kprime = np.min(kprime_distances)
    if min_dist_kprime < tolerance:
        valley_points['Kprime'] = kprime_idx_actual
        valley_points['Kprime_coord'] = kpoints[kprime_idx_actual]

    return valley_points


def analyze_valley_band_structure(
    energies: np.ndarray,
    projections: np.ndarray,
    valley_indices: Dict[str, int],
    efermi: float,
    n_bands_near_fermi: int = 4
) -> Dict:
    """
    Analyze band structure at valley points

    Args:
        energies: Band energies [nkpts, nbands]
        projections: Orbital projections [nkpts, nbands, norbitals]
        valley_indices: Dictionary with valley k-point indices
        efermi: Fermi energy
        n_bands_near_fermi: Number of bands to analyze near Fermi level

    Returns:
        Dictionary with valley analysis results
    """
    results = {}

    energies_shifted = energies - efermi

    for valley_name in ['Gamma', 'K', 'Kprime']:
        if valley_name not in valley_indices:
            continue

        k_idx = valley_indices[valley_name]

        # Get band energies at this k-point
        valley_energies = energies_shifted[k_idx, :]

        # Find bands near Fermi level
        occupied = valley_energies < 0
        unoccupied = valley_energies > 0

        if np.any(occupied):
            vbm_band = np.where(occupied)[0][-1]
            vbm_energy = valley_energies[vbm_band]
        else:
            vbm_band = None
            vbm_energy = None

        if np.any(unoccupied):
            cbm_band = np.where(unoccupied)[0][0]
            cbm_energy = valley_energies[cbm_band]
        else:
            cbm_band = None
            cbm_energy = None

        # Store results
        results[valley_name] = {
            'k_index': k_idx,
            'k_coord': valley_indices.get(f'{valley_name}_coord', None),
            'vbm_band': vbm_band,
            'vbm_energy': vbm_energy,
            'cbm_band': cbm_band,
            'cbm_energy': cbm_energy,
            'all_energies': valley_energies
        }

        # If projections available, analyze orbital character
        if projections is not None and len(projections.shape) == 3:
            # Sum over all atoms
            valley_proj = projections[k_idx, :, :].sum(axis=1)  # [nbands, norbitals] -> [nbands]
            results[valley_name]['orbital_character'] = valley_proj

    return results


def compute_valley_splitting(valley_results: Dict) -> Dict:
    """
    Compute valley splitting between K and K' points

    Args:
        valley_results: Results from analyze_valley_band_structure

    Returns:
        Dictionary with splitting information
    """
    splitting = {}

    if 'K' in valley_results and 'Kprime' in valley_results:
        k_vbm = valley_results['K'].get('vbm_energy')
        kp_vbm = valley_results['Kprime'].get('vbm_energy')

        if k_vbm is not None and kp_vbm is not None:
            splitting['vbm_splitting'] = abs(k_vbm - kp_vbm)
            splitting['vbm_K'] = k_vbm
            splitting['vbm_Kprime'] = kp_vbm

        k_cbm = valley_results['K'].get('cbm_energy')
        kp_cbm = valley_results['Kprime'].get('cbm_energy')

        if k_cbm is not None and kp_cbm is not None:
            splitting['cbm_splitting'] = abs(k_cbm - kp_cbm)
            splitting['cbm_K'] = k_cbm
            splitting['cbm_Kprime'] = kp_cbm

    return splitting


def plot_valley_band_energies(
    valley_results: Dict,
    output: str = 'valley_energies.png',
    title: str = 'Valley Band Energies'
):
    """
    Plot band energies at different valley points

    Args:
        valley_results: Results from analyze_valley_band_structure
        output: Output filename
        title: Plot title
    """
    fig, ax = plt.subplots(figsize=(10, 8))

    valley_names = ['Gamma', 'K', 'Kprime']
    x_positions = [0, 1, 2]
    labels = ['Γ', 'K', "K'"]

    colors = {'Gamma': 'blue', 'K': 'red', 'Kprime': 'green'}

    # Plot band energies
    for valley_name, x_pos, label in zip(valley_names, x_positions, labels):
        if valley_name not in valley_results:
            continue

        energies = valley_results[valley_name]['all_energies']

        # Plot all bands
        for energy in energies:
            ax.plot([x_pos - 0.1, x_pos + 0.1], [energy, energy],
                   color=colors[valley_name], linewidth=2, alpha=0.6)

        # Highlight VBM and CBM
        vbm_energy = valley_results[valley_name].get('vbm_energy')
        cbm_energy = valley_results[valley_name].get('cbm_energy')

        if vbm_energy is not None:
            ax.plot([x_pos - 0.15, x_pos + 0.15], [vbm_energy, vbm_energy],
                   color=colors[valley_name], linewidth=4, label=f'{label} VBM' if label == 'Γ' else '')

        if cbm_energy is not None:
            ax.plot([x_pos - 0.15, x_pos + 0.15], [cbm_energy, cbm_energy],
                   color=colors[valley_name], linewidth=4, linestyle='--',
                   label=f'{label} CBM' if label == 'Γ' else '')

    # Fermi level
    ax.axhline(0, color='k', linestyle='--', linewidth=1.5, alpha=0.7, label='$E_F$')

    # Formatting
    ax.set_xticks(x_positions[:len([v for v in valley_names if v in valley_results])])
    ax.set_xticklabels([labels[i] for i, v in enumerate(valley_names) if v in valley_results],
                       fontsize=16)
    ax.set_ylabel('Energy - $E_F$ (eV)', fontsize=14, fontweight='bold')
    ax.set_title(title, fontsize=16, fontweight='bold')
    ax.set_ylim(-3, 3)
    ax.grid(True, alpha=0.3, axis='y')
    ax.tick_params(labelsize=12)
    ax.legend(fontsize=11, loc='upper right')

    plt.tight_layout()
    plt.savefig(output, dpi=300, bbox_inches='tight')
    print(f"Valley band energies plot saved to: {output}")
    plt.close()


def analyze_orbital_character_at_valleys(
    data: pp.ProcarData,
    valley_indices: Dict[str, int],
    orbital_names: List[str] = None
) -> Dict:
    """
    Analyze orbital character at K and K' valleys

    Focus on d-orbitals for TMDC (d_z2, d_xy, d_x2-y2)

    Args:
        data: ProcarData object
        valley_indices: Valley k-point indices
        orbital_names: Specific orbitals to analyze

    Returns:
        Dictionary with orbital character at each valley
    """
    if orbital_names is None:
        # Focus on d-orbitals for TMDC
        orbital_names = ['dz2', 'dxy', 'dx2']

    orbital_analysis = {}

    energies = data.energies[:, :, 0] if data.is_spin_polarized else data.energies

    for valley_name in ['K', 'Kprime']:
        if valley_name not in valley_indices:
            continue

        k_idx = valley_indices[valley_name]

        # Get VBM and CBM band indices
        valley_energies = energies[k_idx, :]
        occupied = valley_energies < 0
        unoccupied = valley_energies > 0

        if np.any(occupied):
            vbm_band = np.where(occupied)[0][-1]
        else:
            vbm_band = None

        if np.any(unoccupied):
            cbm_band = np.where(unoccupied)[0][0]
        else:
            cbm_band = None

        orbital_analysis[valley_name] = {}

        # Analyze orbital character for VBM and CBM
        for band_type, band_idx in [('VBM', vbm_band), ('CBM', cbm_band)]:
            if band_idx is None:
                continue

            orbital_analysis[valley_name][band_type] = {}

            for orb_name in orbital_names:
                if orb_name in data.orbital_names:
                    orb_idx = data.orbital_names.index(orb_name)

                    if data.is_spin_polarized:
                        # Sum over atoms and spins
                        proj = data.projections[k_idx, band_idx, :, orb_idx, :].sum()
                    else:
                        # Sum over atoms
                        proj = data.projections[k_idx, band_idx, :, orb_idx].sum()

                    orbital_analysis[valley_name][band_type][orb_name] = proj

    return orbital_analysis


def plot_brillouin_zone_with_valleys(
    valley_indices: Dict[str, int],
    output: str = 'brillouin_zone_valleys.png'
):
    """
    Plot hexagonal Brillouin zone with marked valley points

    Args:
        valley_indices: Valley k-point indices and coordinates
        output: Output filename
    """
    fig, ax = plt.subplots(figsize=(8, 8))

    # Draw hexagonal Brillouin zone
    # For hexagonal lattice in 2D
    theta = np.linspace(0, 2*np.pi, 7)
    bz_radius = 2 * np.pi / (3**0.5)  # Simplified

    hex_x = bz_radius * np.cos(theta + np.pi/6)
    hex_y = bz_radius * np.sin(theta + np.pi/6)

    ax.plot(hex_x, hex_y, 'k-', linewidth=2)
    ax.fill(hex_x, hex_y, color='lightgray', alpha=0.3)

    # Mark high-symmetry points
    points_to_plot = {
        'Gamma': (0, 0, 'Γ', 'blue'),
        'K': (2*np.pi/(3*3**0.5), 2*np.pi/3, 'K', 'red'),
        'Kprime': (-2*np.pi/(3*3**0.5), -2*np.pi/3, "K'", 'green')
    }

    for valley_name, (x, y, label, color) in points_to_plot.items():
        if valley_name in valley_indices:
            ax.plot(x, y, 'o', markersize=15, color=color, zorder=10)
            ax.text(x, y - 0.3, label, fontsize=16, fontweight='bold',
                   ha='center', va='top', color=color)

    # Labels
    ax.set_xlabel('$k_x$ (2π/a)', fontsize=14, fontweight='bold')
    ax.set_ylabel('$k_y$ (2π/a)', fontsize=14, fontweight='bold')
    ax.set_title('Hexagonal Brillouin Zone - Valley Points', fontsize=16, fontweight='bold')
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.3)
    ax.tick_params(labelsize=12)

    plt.tight_layout()
    plt.savefig(output, dpi=300, bbox_inches='tight')
    print(f"Brillouin zone plot saved to: {output}")
    plt.close()


def generate_valley_analysis_report(
    valley_results: Dict,
    splitting: Dict,
    orbital_analysis: Dict,
    output_file: str = 'valley_analysis_report.txt'
):
    """
    Generate comprehensive valley analysis report

    Args:
        valley_results: Valley band structure results
        splitting: Valley splitting information
        orbital_analysis: Orbital character at valleys
        output_file: Output text file
    """
    report = []

    report.append("=" * 80)
    report.append("VALLEYTRONICS ANALYSIS REPORT")
    report.append("TMDC Janus Heterostructure - MoSSe/WSSe")
    report.append("=" * 80)

    # Valley energies
    report.append("\n" + "=" * 80)
    report.append("1. VALLEY BAND ENERGIES")
    report.append("=" * 80)

    for valley_name in ['Gamma', 'K', 'Kprime']:
        if valley_name not in valley_results:
            continue

        label = {'Gamma': 'Γ', 'K': 'K', 'Kprime': "K'"}[valley_name]
        result = valley_results[valley_name]

        report.append(f"\n{label} point:")
        report.append(f"  k-point index: {result['k_index']}")

        if result.get('k_coord') is not None:
            k_coord = result['k_coord']
            report.append(f"  k-coordinates: ({k_coord[0]:.4f}, {k_coord[1]:.4f}, {k_coord[2]:.4f})")

        if result['vbm_energy'] is not None:
            report.append(f"  VBM energy: {result['vbm_energy']:.4f} eV")
            report.append(f"  VBM band index: {result['vbm_band']}")

        if result['cbm_energy'] is not None:
            report.append(f"  CBM energy: {result['cbm_energy']:.4f} eV")
            report.append(f"  CBM band index: {result['cbm_band']}")

        if result['vbm_energy'] is not None and result['cbm_energy'] is not None:
            gap = result['cbm_energy'] - result['vbm_energy']
            report.append(f"  Direct gap at {label}: {gap:.4f} eV")

    # Valley splitting
    if splitting:
        report.append("\n" + "=" * 80)
        report.append("2. VALLEY SPLITTING (K vs K')")
        report.append("=" * 80)

        if 'vbm_splitting' in splitting:
            report.append(f"\nValence Band (VBM):")
            report.append(f"  Energy at K:  {splitting['vbm_K']:.4f} eV")
            report.append(f"  Energy at K': {splitting['vbm_Kprime']:.4f} eV")
            report.append(f"  Valley splitting: {splitting['vbm_splitting']:.4f} eV")

            if splitting['vbm_splitting'] < 0.01:
                report.append(f"  → Valleys are nearly degenerate (good for valleytronics)")
            else:
                report.append(f"  → Significant valley splitting detected")

        if 'cbm_splitting' in splitting:
            report.append(f"\nConduction Band (CBM):")
            report.append(f"  Energy at K:  {splitting['cbm_K']:.4f} eV")
            report.append(f"  Energy at K': {splitting['cbm_Kprime']:.4f} eV")
            report.append(f"  Valley splitting: {splitting['cbm_splitting']:.4f} eV")

            if splitting['cbm_splitting'] < 0.01:
                report.append(f"  → Valleys are nearly degenerate (good for valleytronics)")
            else:
                report.append(f"  → Significant valley splitting detected")

    # Orbital character
    if orbital_analysis:
        report.append("\n" + "=" * 80)
        report.append("3. ORBITAL CHARACTER AT VALLEY POINTS")
        report.append("=" * 80)

        for valley_name in ['K', 'Kprime']:
            if valley_name not in orbital_analysis:
                continue

            label = {'K': 'K', 'Kprime': "K'"}[valley_name]
            report.append(f"\n{label} valley:")

            for band_type in ['VBM', 'CBM']:
                if band_type not in orbital_analysis[valley_name]:
                    continue

                report.append(f"\n  {band_type} orbital composition:")
                orb_data = orbital_analysis[valley_name][band_type]

                total = sum(orb_data.values())
                for orb_name, proj in sorted(orb_data.items(), key=lambda x: -x[1]):
                    percentage = 100 * proj / total if total > 0 else 0
                    report.append(f"    {orb_name:8s}: {proj:.4f} ({percentage:.1f}%)")

    # Valleytronics applications
    report.append("\n" + "=" * 80)
    report.append("4. VALLEYTRONICS POTENTIAL")
    report.append("=" * 80)

    if splitting:
        vbm_split = splitting.get('vbm_splitting', float('inf'))
        cbm_split = splitting.get('cbm_splitting', float('inf'))

        report.append(f"\nValley degeneracy assessment:")

        if vbm_split < 0.01 and cbm_split < 0.01:
            report.append("  ✓ Excellent valley degeneracy")
            report.append("  → Strong potential for valley-based information storage")
            report.append("  → Suitable for valley Hall effect devices")
        elif vbm_split < 0.05 or cbm_split < 0.05:
            report.append("  ~ Good valley characteristics")
            report.append("  → Moderate valley splitting may enable valley-selective excitation")
        else:
            report.append("  ✗ Significant valley splitting")
            report.append("  → May reduce valleytronic performance")
            report.append("  → Consider SOC calculations for spin-valley coupling")

    report.append("\n" + "=" * 80)

    # Write to file
    report_text = "\n".join(report)
    with open(output_file, 'w') as f:
        f.write(report_text)

    print(f"\nValley analysis report saved to: {output_file}")
    print(report_text)


def main():
    """Main function"""
    import argparse

    parser = argparse.ArgumentParser(description='Valleytronics analysis for TMDC heterostructures')
    parser.add_argument('--procar', default='PROCAR', help='PROCAR file')
    parser.add_argument('--doscar', default='DOSCAR', help='DOSCAR file')
    parser.add_argument('--output-dir', default='.', help='Output directory')

    args = parser.parse_args()

    # Parse PROCAR
    print(f"Parsing PROCAR: {args.procar}")
    parser_obj = pp.ProcarParser(args.procar)
    data = parser_obj.parse()

    # Read Fermi energy
    from plot_fatband import read_fermi_energy
    efermi = read_fermi_energy(args.doscar)
    print(f"Fermi energy: {efermi:.4f} eV")

    # Identify valley points
    print("\nIdentifying valley points...")
    valley_indices = identify_valley_points(data.kpoints)

    print("Valley points found:")
    for valley_name, idx in valley_indices.items():
        if not valley_name.endswith('_coord'):
            print(f"  {valley_name}: k-point index {idx}")

    # Analyze valley band structure
    energies = data.energies[:, :, 0] if data.is_spin_polarized else data.energies

    print("\nAnalyzing valley band structure...")
    valley_results = analyze_valley_band_structure(
        energies, None, valley_indices, efermi
    )

    # Compute valley splitting
    splitting = compute_valley_splitting(valley_results)

    # Analyze orbital character
    print("\nAnalyzing orbital character at valleys...")
    orbital_analysis = analyze_orbital_character_at_valleys(
        data, valley_indices
    )

    # Generate plots
    print("\nGenerating plots...")

    plot_valley_band_energies(
        valley_results,
        output=os.path.join(args.output_dir, 'valley_energies.png')
    )

    plot_brillouin_zone_with_valleys(
        valley_indices,
        output=os.path.join(args.output_dir, 'brillouin_zone_valleys.png')
    )

    # Generate report
    generate_valley_analysis_report(
        valley_results, splitting, orbital_analysis,
        output_file=os.path.join(args.output_dir, 'valley_analysis_report.txt')
    )

    print("\nValley analysis complete!")


if __name__ == '__main__':
    main()
