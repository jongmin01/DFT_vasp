#!/usr/bin/env python3
"""
plot_fatband.py - Plot orbital-resolved (fat) band structure

Visualizes band structure with orbital character through:
- Line thickness representing orbital contribution
- Color-coded orbital projections
- Multiple orbital overlays

For TMDC Janus heterostructures
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize
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
        print("Warning: Could not read Fermi energy from DOSCAR, using 0.0 eV")
        return 0.0


def find_high_symmetry_points(kpoints: np.ndarray, threshold: float = 0.01) -> Tuple[List[int], List[str]]:
    """
    Automatically detect high-symmetry points in k-path

    Args:
        kpoints: k-point coordinates
        threshold: Threshold for detecting kink in k-path

    Returns:
        indices: List of k-point indices
        labels: List of labels (auto-generated)
    """
    kdist = get_kpath_distance(kpoints)
    high_sym_indices = [0]  # Always include first point

    # Detect kinks in k-path
    for i in range(1, len(kpoints) - 1):
        # Check if direction changes
        v1 = kpoints[i] - kpoints[i-1]
        v2 = kpoints[i+1] - kpoints[i]

        # Normalize
        v1_norm = np.linalg.norm(v1)
        v2_norm = np.linalg.norm(v2)

        if v1_norm > 1e-8 and v2_norm > 1e-8:
            v1 = v1 / v1_norm
            v2 = v2 / v2_norm

            # Check angle
            cos_angle = np.dot(v1, v2)
            if cos_angle < 0.99:  # Significant direction change
                high_sym_indices.append(i)

    high_sym_indices.append(len(kpoints) - 1)  # Always include last point

    # Generate labels (for hexagonal system, common for TMDC)
    # Try to identify Gamma, M, K points
    labels = []
    for idx in high_sym_indices:
        kpt = kpoints[idx]
        # Check if it's Gamma point
        if np.allclose(kpt, [0, 0, 0], atol=0.01):
            labels.append('Γ')
        else:
            labels.append('')  # Will be set based on position

    # Common TMDC path: Γ-M-K-Γ
    if len(labels) == 4:
        labels = ['Γ', 'M', 'K', 'Γ']
    elif len(labels) == 3:
        labels = ['Γ', 'K', 'Γ']

    return high_sym_indices, labels


def plot_fatband_single_orbital(
    kpoints: np.ndarray,
    energies: np.ndarray,
    projections: np.ndarray,
    orbital_name: str,
    efermi: float = 0.0,
    high_sym_points: Optional[List[int]] = None,
    high_sym_labels: Optional[List[str]] = None,
    title: str = 'Fat Band Structure',
    output: str = 'fatband.png',
    energy_range: Tuple[float, float] = (-3, 3),
    projection_scale: float = 20.0
):
    """
    Plot fat band structure for a single orbital

    Args:
        kpoints: k-point coordinates [nkpts, 3]
        energies: Band energies [nkpts, nbands]
        projections: Orbital projections [nkpts, nbands]
        orbital_name: Name of orbital being plotted
        efermi: Fermi energy
        high_sym_points: High-symmetry point indices
        high_sym_labels: High-symmetry point labels
        title: Plot title
        output: Output filename
        energy_range: Energy window (relative to Fermi level)
        projection_scale: Scaling factor for projection visualization
    """
    kdist = get_kpath_distance(kpoints)
    energies_shifted = energies - efermi

    # Auto-detect high-symmetry points if not provided
    if high_sym_points is None:
        high_sym_points, high_sym_labels = find_high_symmetry_points(kpoints)

    # Create figure
    fig, ax = plt.subplots(figsize=(10, 8))

    # Plot bands with width proportional to projection
    for iband in range(energies.shape[1]):
        # Get energies and projections for this band
        band_energies = energies_shifted[:, iband]
        band_proj = projections[:, iband]

        # Filter by energy range
        mask = (band_energies >= energy_range[0]) & (band_energies <= energy_range[1])

        if not np.any(mask):
            continue

        # Plot with varying linewidth
        for i in range(len(kdist) - 1):
            if not (mask[i] or mask[i+1]):
                continue

            # Line width based on projection
            lw = 0.5 + projection_scale * (band_proj[i] + band_proj[i+1]) / 2

            ax.plot(
                kdist[i:i+2],
                band_energies[i:i+2],
                color='blue',
                linewidth=lw,
                alpha=0.6,
                solid_capstyle='round'
            )

    # Fermi level
    ax.axhline(0, color='k', linestyle='--', linewidth=1.5, alpha=0.7, label='$E_F$')

    # High-symmetry points
    if high_sym_points is not None and high_sym_labels is not None:
        for kpt in high_sym_points:
            ax.axvline(kdist[kpt], color='gray', linestyle='-', linewidth=0.8, alpha=0.5)

        ax.set_xticks([kdist[kpt] for kpt in high_sym_points])
        ax.set_xticklabels(high_sym_labels, fontsize=14)

    # Labels
    ax.set_ylabel('Energy - $E_F$ (eV)', fontsize=14, fontweight='bold')
    ax.set_title(f'{title}\nOrbital: {orbital_name}', fontsize=16, fontweight='bold')
    ax.set_xlim(kdist[0], kdist[-1])
    ax.set_ylim(energy_range)
    ax.legend(fontsize=12, loc='upper right')
    ax.grid(True, alpha=0.3, linestyle='--', axis='y')
    ax.tick_params(labelsize=12)

    plt.tight_layout()
    plt.savefig(output, dpi=300, bbox_inches='tight')
    print(f"Fat band plot saved to: {output}")
    plt.close()


def plot_fatband_multiorbital(
    kpoints: np.ndarray,
    energies: np.ndarray,
    projections_dict: Dict[str, np.ndarray],
    efermi: float = 0.0,
    high_sym_points: Optional[List[int]] = None,
    high_sym_labels: Optional[List[str]] = None,
    title: str = 'Fat Band Structure',
    output: str = 'fatband_multi.png',
    energy_range: Tuple[float, float] = (-3, 3),
    cmap: str = 'rainbow'
):
    """
    Plot fat band structure with multiple orbitals (color-coded)

    Args:
        kpoints: k-point coordinates [nkpts, 3]
        energies: Band energies [nkpts, nbands]
        projections_dict: Dictionary of orbital projections {orbital_name: projection_array}
        efermi: Fermi energy
        high_sym_points: High-symmetry point indices
        high_sym_labels: High-symmetry point labels
        title: Plot title
        output: Output filename
        energy_range: Energy window (relative to Fermi level)
        cmap: Colormap for different orbitals
    """
    kdist = get_kpath_distance(kpoints)
    energies_shifted = energies - efermi

    # Auto-detect high-symmetry points if not provided
    if high_sym_points is None:
        high_sym_points, high_sym_labels = find_high_symmetry_points(kpoints)

    # Create figure
    fig, ax = plt.subplots(figsize=(10, 8))

    # Assign colors to orbitals
    orbital_names = list(projections_dict.keys())
    colors = plt.cm.get_cmap(cmap)(np.linspace(0, 1, len(orbital_names)))
    orbital_colors = {name: colors[i] for i, name in enumerate(orbital_names)}

    # Plot bands for each orbital
    for orbital_name, color in orbital_colors.items():
        projections = projections_dict[orbital_name]

        for iband in range(energies.shape[1]):
            band_energies = energies_shifted[:, iband]
            band_proj = projections[:, iband]

            # Filter by energy range
            mask = (band_energies >= energy_range[0]) & (band_energies <= energy_range[1])

            if not np.any(mask):
                continue

            # Create scatter plot with size proportional to projection
            scatter_size = 50 * band_proj[mask]

            ax.scatter(
                kdist[mask],
                band_energies[mask],
                s=scatter_size,
                color=color,
                alpha=0.6,
                edgecolors='none',
                label=orbital_name if iband == 0 else ''
            )

    # Fermi level
    ax.axhline(0, color='k', linestyle='--', linewidth=1.5, alpha=0.7, label='$E_F$')

    # High-symmetry points
    if high_sym_points is not None and high_sym_labels is not None:
        for kpt in high_sym_points:
            ax.axvline(kdist[kpt], color='gray', linestyle='-', linewidth=0.8, alpha=0.5)

        ax.set_xticks([kdist[kpt] for kpt in high_sym_points])
        ax.set_xticklabels(high_sym_labels, fontsize=14)

    # Labels
    ax.set_ylabel('Energy - $E_F$ (eV)', fontsize=14, fontweight='bold')
    ax.set_title(title, fontsize=16, fontweight='bold')
    ax.set_xlim(kdist[0], kdist[-1])
    ax.set_ylim(energy_range)

    # Legend (remove duplicates)
    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ax.legend(by_label.values(), by_label.keys(), fontsize=11, loc='upper right',
              framealpha=0.9, ncol=2)

    ax.grid(True, alpha=0.3, linestyle='--', axis='y')
    ax.tick_params(labelsize=12)

    plt.tight_layout()
    plt.savefig(output, dpi=300, bbox_inches='tight')
    print(f"Multi-orbital fat band plot saved to: {output}")
    plt.close()


def plot_fatband_colormap(
    kpoints: np.ndarray,
    energies: np.ndarray,
    projections: np.ndarray,
    projection_name: str,
    efermi: float = 0.0,
    high_sym_points: Optional[List[int]] = None,
    high_sym_labels: Optional[List[str]] = None,
    title: str = 'Fat Band Structure',
    output: str = 'fatband_colormap.png',
    energy_range: Tuple[float, float] = (-3, 3),
    cmap: str = 'hot'
):
    """
    Plot fat band structure using colormap for projection intensity

    Args:
        kpoints: k-point coordinates [nkpts, 3]
        energies: Band energies [nkpts, nbands]
        projections: Projections [nkpts, nbands]
        projection_name: Name of projection (e.g., "Mo d-orbitals", "MoSSe layer")
        efermi: Fermi energy
        high_sym_points: High-symmetry point indices
        high_sym_labels: High-symmetry point labels
        title: Plot title
        output: Output filename
        energy_range: Energy window (relative to Fermi level)
        cmap: Colormap name
    """
    kdist = get_kpath_distance(kpoints)
    energies_shifted = energies - efermi

    # Auto-detect high-symmetry points if not provided
    if high_sym_points is None:
        high_sym_points, high_sym_labels = find_high_symmetry_points(kpoints)

    # Create figure
    fig, ax = plt.subplots(figsize=(10, 8))

    # Normalize projection for colormap
    norm = Normalize(vmin=0, vmax=projections.max())
    colormap = plt.cm.get_cmap(cmap)

    # Plot bands with color based on projection
    for iband in range(energies.shape[1]):
        band_energies = energies_shifted[:, iband]
        band_proj = projections[:, iband]

        # Filter by energy range
        mask = (band_energies >= energy_range[0]) & (band_energies <= energy_range[1])

        if not np.any(mask):
            continue

        # Create line collection with varying colors
        points = np.array([kdist, band_energies]).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)

        lc = LineCollection(segments, cmap=cmap, norm=norm, linewidths=2, alpha=0.8)
        lc.set_array(band_proj)
        ax.add_collection(lc)

    # Fermi level
    ax.axhline(0, color='k', linestyle='--', linewidth=1.5, alpha=0.7, label='$E_F$')

    # High-symmetry points
    if high_sym_points is not None and high_sym_labels is not None:
        for kpt in high_sym_points:
            ax.axvline(kdist[kpt], color='gray', linestyle='-', linewidth=0.8, alpha=0.5)

        ax.set_xticks([kdist[kpt] for kpt in high_sym_points])
        ax.set_xticklabels(high_sym_labels, fontsize=14)

    # Colorbar
    sm = plt.cm.ScalarMappable(cmap=colormap, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, pad=0.02)
    cbar.set_label(f'{projection_name} projection', fontsize=12, fontweight='bold')
    cbar.ax.tick_params(labelsize=11)

    # Labels
    ax.set_ylabel('Energy - $E_F$ (eV)', fontsize=14, fontweight='bold')
    ax.set_title(title, fontsize=16, fontweight='bold')
    ax.set_xlim(kdist[0], kdist[-1])
    ax.set_ylim(energy_range)
    ax.legend(fontsize=12, loc='upper left')
    ax.grid(True, alpha=0.3, linestyle='--', axis='y')
    ax.tick_params(labelsize=12)

    plt.tight_layout()
    plt.savefig(output, dpi=300, bbox_inches='tight')
    print(f"Colormap fat band plot saved to: {output}")
    plt.close()


def main():
    """Main function for standalone usage"""
    import argparse

    parser = argparse.ArgumentParser(description='Plot orbital-resolved (fat) band structure')
    parser.add_argument('--procar', default='PROCAR', help='PROCAR file path')
    parser.add_argument('--doscar', default='DOSCAR', help='DOSCAR file path (for Fermi energy)')
    parser.add_argument('--orbital', type=str, help='Orbital name for single-orbital plot')
    parser.add_argument('--orbitals', nargs='+', help='List of orbitals for multi-orbital plot')
    parser.add_argument('--output', default='fatband.png', help='Output filename')
    parser.add_argument('--title', default='Fat Band Structure', help='Plot title')
    parser.add_argument('--emin', type=float, default=-3.0, help='Minimum energy (eV)')
    parser.add_argument('--emax', type=float, default=3.0, help='Maximum energy (eV)')
    parser.add_argument('--style', choices=['single', 'multi', 'colormap'], default='single',
                       help='Plot style')

    args = parser.parse_args()

    # Parse PROCAR
    print(f"Parsing PROCAR: {args.procar}")
    parser_obj = pp.ProcarParser(args.procar)
    data = parser_obj.parse()

    print(f"  K-points: {data.nkpts}")
    print(f"  Bands: {data.nbands}")
    print(f"  Ions: {data.nions}")
    print(f"  Orbitals: {data.orbital_names}")

    # Read Fermi energy
    efermi = read_fermi_energy(args.doscar)
    print(f"  Fermi energy: {efermi:.4f} eV")

    # Get orbital projections
    if args.style == 'single':
        if args.orbital is None:
            print("Error: --orbital required for single-orbital plot")
            sys.exit(1)

        proj = parser_obj.get_orbital_projection(args.orbital)
        if data.is_spin_polarized:
            proj = proj[:, :, :, 0]  # Use spin-up
        proj_sum = proj.sum(axis=2)  # Sum over atoms

        plot_fatband_single_orbital(
            data.kpoints, data.energies[:, :, 0] if data.is_spin_polarized else data.energies,
            proj_sum, args.orbital, efermi,
            title=args.title, output=args.output,
            energy_range=(args.emin, args.emax)
        )

    elif args.style == 'multi':
        if args.orbitals is None:
            print("Error: --orbitals required for multi-orbital plot")
            sys.exit(1)

        proj_dict = {}
        for orb in args.orbitals:
            proj = parser_obj.get_orbital_projection(orb)
            if data.is_spin_polarized:
                proj = proj[:, :, :, 0]
            proj_dict[orb] = proj.sum(axis=2)

        plot_fatband_multiorbital(
            data.kpoints, data.energies[:, :, 0] if data.is_spin_polarized else data.energies,
            proj_dict, efermi,
            title=args.title, output=args.output,
            energy_range=(args.emin, args.emax)
        )

    print("\nDone!")


if __name__ == '__main__':
    main()
