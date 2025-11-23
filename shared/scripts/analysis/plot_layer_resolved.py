#!/usr/bin/env python3
"""
plot_layer_resolved.py - Layer-resolved band structure for heterostructures

Analyzes and visualizes layer contributions to band structure:
- Layer-projected band structure
- Interlayer coupling analysis
- Band alignment (Type-I/Type-II)
- Charge transfer identification

Specifically designed for TMDC Janus heterostructures (MoSSe/WSSe)
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import sys
import os
from typing import List, Dict, Optional, Tuple

# Import local modules
import procar_parser as pp


def get_kpath_distance(kpoints: np.ndarray) -> np.ndarray:
    """Calculate cumulative distance along k-path"""
    distances = [0.0]
    for i in range(1, len(kpoints)):
        dk = np.linalg.norm(kpoints[i] - kpoints[i-1])
        distances.append(distances[-1] + dk)
    return np.array(distances)


def read_fermi_energy(doscar_file: str = 'DOSCAR') -> float:
    """Read Fermi energy from DOSCAR"""
    try:
        with open(doscar_file, 'r') as f:
            lines = f.readlines()
        header = lines[5].split()
        efermi = float(header[3])
        return efermi
    except:
        return 0.0


def identify_tmdc_layers(poscar_file: str) -> Dict[str, List[int]]:
    """
    Identify TMDC layers from POSCAR file

    For MoSSe/WSSe heterostructure, atoms are typically ordered as:
    Mo S Se W S Se (or similar)

    Returns:
        Dictionary mapping layer names to atom indices
    """
    with open(poscar_file, 'r') as f:
        lines = f.readlines()

    # Line 5: atom types
    # Line 6: number of each atom type
    atom_types = lines[5].split()
    atom_counts = [int(x) for x in lines[6].split()]

    print(f"POSCAR atom types: {atom_types}")
    print(f"POSCAR atom counts: {atom_counts}")

    # Create atom index mapping
    atom_indices = {}
    current_idx = 0

    for atom_type, count in zip(atom_types, atom_counts):
        indices = list(range(current_idx, current_idx + count))
        atom_indices[atom_type] = indices
        current_idx += count

    # Identify layers based on atom types
    # Typical ordering for MoSSe/WSSe: Mo, S(top), Se(bottom), W, S(top), Se(bottom)
    # Or could be grouped differently

    layers = {}

    # Try to identify MoSSe and WSSe layers
    if 'Mo' in atom_indices and 'W' in atom_indices:
        # Find which S/Se belong to which layer
        mo_indices = atom_indices['Mo']
        w_indices = atom_indices['W']
        s_indices = atom_indices.get('S', [])
        se_indices = atom_indices.get('Se', [])

        # Assume first half of S/Se belong to MoSSe, second half to WSSe
        n_s = len(s_indices)
        n_se = len(se_indices)

        # MoSSe layer: Mo + half S + half Se
        mosse_s = s_indices[:n_s//2] if n_s > 0 else []
        mosse_se = se_indices[:n_se//2] if n_se > 0 else []
        layers['MoSSe'] = mo_indices + mosse_s + mosse_se

        # WSSe layer: W + half S + half Se
        wsse_s = s_indices[n_s//2:] if n_s > 0 else []
        wsse_se = se_indices[n_se//2:] if n_se > 0 else []
        layers['WSSe'] = w_indices + wsse_s + wsse_se

        # Sublayers
        layers['Mo'] = mo_indices
        layers['W'] = w_indices
        layers['S_MoSSe'] = mosse_s
        layers['Se_MoSSe'] = mosse_se
        layers['S_WSSe'] = wsse_s
        layers['Se_WSSe'] = wsse_se

    print(f"\nIdentified layers:")
    for layer_name, indices in layers.items():
        print(f"  {layer_name}: atoms {indices}")

    return layers


def find_vbm_cbm(energies: np.ndarray, efermi: float) -> Tuple:
    """
    Find VBM and CBM positions

    Returns:
        vbm_k, vbm_band, vbm_energy, cbm_k, cbm_band, cbm_energy, gap
    """
    energies_shifted = energies - efermi

    # Find VBM (highest occupied state)
    occupied = energies_shifted < 0
    if not np.any(occupied):
        return None, None, None, None, None, None, None

    vbm_indices = np.where(occupied)
    vbm_flat = energies_shifted[occupied].argmax()
    vbm_k = vbm_indices[0][vbm_flat]
    vbm_band = vbm_indices[1][vbm_flat]
    vbm_energy = energies_shifted[vbm_k, vbm_band]

    # Find CBM (lowest unoccupied state)
    unoccupied = energies_shifted > 0
    if not np.any(unoccupied):
        return vbm_k, vbm_band, vbm_energy, None, None, None, None

    cbm_indices = np.where(unoccupied)
    cbm_flat = energies_shifted[unoccupied].argmin()
    cbm_k = cbm_indices[0][cbm_flat]
    cbm_band = cbm_indices[1][cbm_flat]
    cbm_energy = energies_shifted[cbm_k, cbm_band]

    gap = cbm_energy - vbm_energy

    return vbm_k, vbm_band, vbm_energy, cbm_k, cbm_band, cbm_energy, gap


def plot_layer_resolved_bands(
    kpoints: np.ndarray,
    energies: np.ndarray,
    layer_projections: Dict[str, np.ndarray],
    efermi: float = 0.0,
    high_sym_points: Optional[List[int]] = None,
    high_sym_labels: Optional[List[str]] = None,
    title: str = 'Layer-Resolved Band Structure',
    output: str = 'layer_bands.png',
    energy_range: Tuple[float, float] = (-3, 3)
):
    """
    Plot layer-resolved band structure

    Args:
        kpoints: k-point coordinates [nkpts, 3]
        energies: Band energies [nkpts, nbands]
        layer_projections: Dictionary of layer projections {layer_name: projection_array [nkpts, nbands]}
        efermi: Fermi energy
        high_sym_points: High-symmetry point indices
        high_sym_labels: High-symmetry point labels
        title: Plot title
        output: Output filename
        energy_range: Energy window (relative to Fermi level)
    """
    kdist = get_kpath_distance(kpoints)
    energies_shifted = energies - efermi

    # Auto-detect high-symmetry points if not provided
    if high_sym_points is None:
        from plot_fatband import find_high_symmetry_points
        high_sym_points, high_sym_labels = find_high_symmetry_points(kpoints)

    # Create figure with multiple subplots
    n_layers = len(layer_projections)
    fig, axes = plt.subplots(1, n_layers + 1, figsize=(5 * (n_layers + 1), 8),
                            sharey=True, gridspec_kw={'width_ratios': [1] * n_layers + [1]})

    if n_layers == 1:
        axes = [axes]

    layer_names = list(layer_projections.keys())
    colors = plt.cm.get_cmap('tab10')(np.linspace(0, 1, n_layers))

    # Plot each layer in separate subplot
    for idx, (layer_name, color) in enumerate(zip(layer_names, colors)):
        ax = axes[idx]
        projections = layer_projections[layer_name]

        # Plot bands with color intensity based on layer projection
        for iband in range(energies.shape[1]):
            band_energies = energies_shifted[:, iband]
            band_proj = projections[:, iband]

            # Filter by energy range
            mask = (band_energies >= energy_range[0]) & (band_energies <= energy_range[1])

            if not np.any(mask):
                continue

            # Plot with varying linewidth/alpha based on projection
            for i in range(len(kdist) - 1):
                if not (mask[i] or mask[i+1]):
                    continue

                # Alpha based on average projection
                alpha = 0.3 + 0.7 * (band_proj[i] + band_proj[i+1]) / (2 * projections.max())

                ax.plot(
                    kdist[i:i+2],
                    band_energies[i:i+2],
                    color=color,
                    linewidth=1.5,
                    alpha=alpha,
                    solid_capstyle='round'
                )

        # Fermi level
        ax.axhline(0, color='k', linestyle='--', linewidth=1.5, alpha=0.7)

        # High-symmetry points
        if high_sym_points is not None:
            for kpt in high_sym_points:
                ax.axvline(kdist[kpt], color='gray', linestyle='-', linewidth=0.8, alpha=0.5)

            ax.set_xticks([kdist[kpt] for kpt in high_sym_points])
            ax.set_xticklabels(high_sym_labels if high_sym_labels else [''] * len(high_sym_points),
                              fontsize=12)

        ax.set_title(f'{layer_name}', fontsize=14, fontweight='bold', color=color)
        ax.set_xlim(kdist[0], kdist[-1])
        ax.set_ylim(energy_range)
        ax.grid(True, alpha=0.3, linestyle='--', axis='y')
        ax.tick_params(labelsize=11)

        if idx == 0:
            ax.set_ylabel('Energy - $E_F$ (eV)', fontsize=13, fontweight='bold')

    # Combined plot in last subplot
    ax = axes[-1]
    for idx, (layer_name, color) in enumerate(zip(layer_names, colors)):
        projections = layer_projections[layer_name]

        for iband in range(energies.shape[1]):
            band_energies = energies_shifted[:, iband]
            band_proj = projections[:, iband]

            mask = (band_energies >= energy_range[0]) & (band_energies <= energy_range[1])
            if not np.any(mask):
                continue

            # Scatter plot
            scatter_size = 30 * band_proj[mask] / projections.max()
            ax.scatter(
                kdist[mask],
                band_energies[mask],
                s=scatter_size,
                color=color,
                alpha=0.7,
                edgecolors='none',
                label=layer_name if iband == 0 else ''
            )

    ax.axhline(0, color='k', linestyle='--', linewidth=1.5, alpha=0.7)

    if high_sym_points is not None:
        for kpt in high_sym_points:
            ax.axvline(kdist[kpt], color='gray', linestyle='-', linewidth=0.8, alpha=0.5)

        ax.set_xticks([kdist[kpt] for kpt in high_sym_points])
        ax.set_xticklabels(high_sym_labels if high_sym_labels else [''] * len(high_sym_points),
                          fontsize=12)

    ax.set_title('Combined', fontsize=14, fontweight='bold')
    ax.set_xlim(kdist[0], kdist[-1])
    ax.set_ylim(energy_range)
    ax.grid(True, alpha=0.3, linestyle='--', axis='y')
    ax.tick_params(labelsize=11)

    # Legend
    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ax.legend(by_label.values(), by_label.keys(), fontsize=11, loc='upper right',
             framealpha=0.9)

    fig.suptitle(title, fontsize=16, fontweight='bold', y=0.98)
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig(output, dpi=300, bbox_inches='tight')
    print(f"Layer-resolved band structure saved to: {output}")
    plt.close()


def analyze_band_character(
    energies: np.ndarray,
    layer_projections: Dict[str, np.ndarray],
    efermi: float,
    output_file: str = 'band_character_analysis.txt'
):
    """
    Analyze VBM and CBM character

    Args:
        energies: Band energies [nkpts, nbands]
        layer_projections: Layer projections
        efermi: Fermi energy
        output_file: Output text file
    """
    # Find VBM and CBM
    vbm_k, vbm_band, vbm_energy, cbm_k, cbm_band, cbm_energy, gap = find_vbm_cbm(energies, efermi)

    if vbm_k is None or cbm_k is None:
        print("Warning: Could not find VBM/CBM")
        return

    # Analyze layer character
    report = []
    report.append("=" * 70)
    report.append("BAND CHARACTER ANALYSIS")
    report.append("=" * 70)
    report.append(f"\nBand Gap: {gap:.4f} eV")
    report.append(f"Gap Type: {'Direct' if vbm_k == cbm_k else 'Indirect'}")
    report.append(f"\nVBM:")
    report.append(f"  Energy: {vbm_energy:.4f} eV (relative to E_F)")
    report.append(f"  k-point index: {vbm_k}")
    report.append(f"  Band index: {vbm_band}")
    report.append(f"\n  Layer contributions:")

    vbm_total = 0
    for layer_name, projections in layer_projections.items():
        vbm_proj = projections[vbm_k, vbm_band]
        vbm_total += vbm_proj
        report.append(f"    {layer_name:15s}: {vbm_proj:8.4f} ({100*vbm_proj:.2f}%)")

    report.append(f"    {'Total':15s}: {vbm_total:8.4f} ({100*vbm_total:.2f}%)")

    report.append(f"\nCBM:")
    report.append(f"  Energy: {cbm_energy:.4f} eV (relative to E_F)")
    report.append(f"  k-point index: {cbm_k}")
    report.append(f"  Band index: {cbm_band}")
    report.append(f"\n  Layer contributions:")

    cbm_total = 0
    for layer_name, projections in layer_projections.items():
        cbm_proj = projections[cbm_k, cbm_band]
        cbm_total += cbm_proj
        report.append(f"    {layer_name:15s}: {cbm_proj:8.4f} ({100*cbm_proj:.2f}%)")

    report.append(f"    {'Total':15s}: {cbm_total:8.4f} ({100*cbm_total:.2f}%)")

    # Determine band alignment type
    report.append(f"\n" + "=" * 70)
    report.append("BAND ALIGNMENT ANALYSIS")
    report.append("=" * 70)

    # Find dominant layers for VBM and CBM
    vbm_dominant = max(layer_projections.items(), key=lambda x: x[1][vbm_k, vbm_band])
    cbm_dominant = max(layer_projections.items(), key=lambda x: x[1][cbm_k, cbm_band])

    report.append(f"\nVBM dominated by: {vbm_dominant[0]} ({100*vbm_dominant[1][vbm_k, vbm_band]:.1f}%)")
    report.append(f"CBM dominated by: {cbm_dominant[0]} ({100*cbm_dominant[1][cbm_k, cbm_band]:.1f}%)")

    # Determine Type-I or Type-II
    if vbm_dominant[0] == cbm_dominant[0]:
        alignment_type = "Type-I (straddling)"
        report.append(f"\nBand Alignment: {alignment_type}")
        report.append("  → VBM and CBM localized on the same layer")
        report.append("  → Charge carriers confined to one layer")
    else:
        alignment_type = "Type-II (staggered)"
        report.append(f"\nBand Alignment: {alignment_type}")
        report.append("  → VBM and CBM localized on different layers")
        report.append("  → Spatial separation of electrons and holes")
        report.append("  → Favorable for charge separation and exciton formation")

    report.append("=" * 70)

    # Write to file
    report_text = "\n".join(report)
    with open(output_file, 'w') as f:
        f.write(report_text)

    print(f"\nBand character analysis saved to: {output_file}")
    print(report_text)


def main():
    """Main function"""
    import argparse

    parser = argparse.ArgumentParser(description='Layer-resolved band structure analysis')
    parser.add_argument('--procar', default='PROCAR', help='PROCAR file')
    parser.add_argument('--poscar', default='POSCAR', help='POSCAR file')
    parser.add_argument('--doscar', default='DOSCAR', help='DOSCAR file')
    parser.add_argument('--layers', nargs='+', help='Layer names to analyze (e.g., MoSSe WSSe)')
    parser.add_argument('--output', default='layer_bands.png', help='Output image file')
    parser.add_argument('--report', default='band_character.txt', help='Analysis report file')
    parser.add_argument('--emin', type=float, default=-3.0, help='Min energy (eV)')
    parser.add_argument('--emax', type=float, default=3.0, help='Max energy (eV)')

    args = parser.parse_args()

    # Parse PROCAR
    print(f"Parsing PROCAR: {args.procar}")
    parser_obj = pp.ProcarParser(args.procar)
    data = parser_obj.parse()

    # Identify layers
    print(f"\nReading structure from: {args.poscar}")
    layer_atoms = identify_tmdc_layers(args.poscar)

    # Get layer projections
    layer_projections = {}
    for layer_name in (args.layers if args.layers else ['MoSSe', 'WSSe']):
        if layer_name in layer_atoms:
            proj = parser_obj.get_atom_projection(layer_atoms[layer_name])
            if data.is_spin_polarized:
                proj = proj[:, :, 0]  # Use spin-up
            layer_projections[layer_name] = proj

    # Read Fermi energy
    efermi = read_fermi_energy(args.doscar)
    print(f"\nFermi energy: {efermi:.4f} eV")

    # Get energies
    energies = data.energies[:, :, 0] if data.is_spin_polarized else data.energies

    # Plot layer-resolved bands
    plot_layer_resolved_bands(
        data.kpoints, energies, layer_projections, efermi,
        title='Layer-Resolved Band Structure - MoSSe/WSSe Heterostructure',
        output=args.output,
        energy_range=(args.emin, args.emax)
    )

    # Analyze band character
    analyze_band_character(energies, layer_projections, efermi, args.report)

    print("\nDone!")


if __name__ == '__main__':
    main()
