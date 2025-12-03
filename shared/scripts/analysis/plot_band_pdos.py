#!/usr/bin/env python3
"""
plot_band_pdos.py - Combined band structure and projected DOS (PDOS) plotting

Creates publication-quality combined plots of:
- Band structure (left panel)
- Projected DOS (right panel)
- Atom/orbital/layer-resolved PDOS

For TMDC Janus heterostructures
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from typing import List, Dict, Optional, Tuple
import sys
import os
import gzip

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
    """Read Fermi energy from DOSCAR (supports .gz compressed files)"""
    try:
        # Check if file is gzip compressed
        if doscar_file.endswith('.gz'):
            with gzip.open(doscar_file, 'rt', encoding='utf-8') as f:
                lines = f.readlines()
        else:
            with open(doscar_file, 'r') as f:
                lines = f.readlines()
        header = lines[5].split()
        efermi = float(header[3])
        return efermi
    except:
        return 0.0


def read_doscar_total(doscar_file: str = 'DOSCAR') -> Tuple[np.ndarray, np.ndarray, float]:
    """
    Read total DOS from DOSCAR (supports .gz compressed files)

    Returns:
        energy: Energy values
        dos: Total DOS
        efermi: Fermi energy
    """
    # Check if file is gzip compressed
    if doscar_file.endswith('.gz'):
        with gzip.open(doscar_file, 'rt', encoding='utf-8') as f:
            lines = f.readlines()
    else:
        with open(doscar_file, 'r') as f:
            lines = f.readlines()

    # Line 6: NEDOS, Fermi energy
    header = lines[5].split()
    nedos = int(header[2])
    efermi = float(header[3])

    # Read DOS data (starts at line 6)
    dos_data = []
    for i in range(6, 6 + nedos):
        dos_data.append([float(x) for x in lines[i].split()])

    dos_data = np.array(dos_data)

    energy = dos_data[:, 0] - efermi  # Shift to Fermi level = 0
    dos = dos_data[:, 1]  # Total DOS (or spin-up for spin-polarized)

    return energy, dos, efermi


def compute_pdos_from_procar(
    data: pp.ProcarData,
    atom_groups: Dict[str, List[int]],
    orbital_groups: Optional[Dict[str, List[str]]] = None
) -> Dict[str, np.ndarray]:
    """
    Compute PDOS from PROCAR projections

    Args:
        data: ProcarData object
        atom_groups: Dictionary of atom groups {group_name: [atom_indices]}
        orbital_groups: Dictionary of orbital groups {group_name: [orbital_names]}

    Returns:
        Dictionary of PDOS arrays
    """
    # This is a simplified version - would need proper DOS calculation
    # For now, we'll use histogram of eigenvalues weighted by projections

    pdos_dict = {}

    # Energy grid
    energies = data.energies[:, :, 0] if data.is_spin_polarized else data.energies
    e_min = energies.min()
    e_max = energies.max()
    n_bins = 1000
    e_grid = np.linspace(e_min, e_max, n_bins)

    # For each atom group
    for group_name, atom_indices in atom_groups.items():
        # Sum projections over atoms and orbitals
        if data.is_spin_polarized:
            proj = data.projections[:, :, atom_indices, :, 0].sum(axis=(2, 3))
        else:
            proj = data.projections[:, :, atom_indices, :].sum(axis=(2, 3))

        # Create histogram
        dos, _ = np.histogram(energies.flatten(), bins=e_grid,
                             weights=proj.flatten(), density=True)

        # Smooth DOS
        from scipy.ndimage import gaussian_filter1d
        dos_smooth = gaussian_filter1d(dos, sigma=2)

        pdos_dict[group_name] = dos_smooth

    e_centers = (e_grid[:-1] + e_grid[1:]) / 2

    return e_centers, pdos_dict


def plot_band_and_pdos(
    kpoints: np.ndarray,
    energies: np.ndarray,
    dos_energy: np.ndarray,
    pdos_dict: Dict[str, np.ndarray],
    efermi: float = 0.0,
    high_sym_points: Optional[List[int]] = None,
    high_sym_labels: Optional[List[str]] = None,
    title: str = 'Band Structure and PDOS',
    output: str = 'band_pdos.png',
    energy_range: Tuple[float, float] = (-3, 3)
):
    """
    Plot combined band structure and PDOS

    Args:
        kpoints: k-point coordinates [nkpts, 3]
        energies: Band energies [nkpts, nbands]
        dos_energy: DOS energy grid
        pdos_dict: Dictionary of PDOS {name: dos_array}
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

    # Create figure with GridSpec
    fig = plt.figure(figsize=(14, 8))
    gs = gridspec.GridSpec(1, 2, width_ratios=[2, 1], wspace=0.05)

    # Left panel: Band structure
    ax_band = fig.add_subplot(gs[0])

    # Plot bands
    for iband in range(energies.shape[1]):
        band_energies = energies_shifted[:, iband]
        mask = (band_energies >= energy_range[0]) & (band_energies <= energy_range[1])

        if np.any(mask):
            ax_band.plot(kdist, band_energies, 'b-', linewidth=1.5, alpha=0.7)

    # Fermi level
    ax_band.axhline(0, color='k', linestyle='--', linewidth=1.5, alpha=0.7, label='$E_F$')

    # High-symmetry points
    if high_sym_points is not None and high_sym_labels is not None:
        for kpt in high_sym_points:
            ax_band.axvline(kdist[kpt], color='gray', linestyle='-', linewidth=0.8, alpha=0.5)

        ax_band.set_xticks([kdist[kpt] for kpt in high_sym_points])
        ax_band.set_xticklabels(high_sym_labels, fontsize=14)

    ax_band.set_ylabel('Energy - $E_F$ (eV)', fontsize=14, fontweight='bold')
    ax_band.set_title('Band Structure', fontsize=14, fontweight='bold')
    ax_band.set_xlim(kdist[0], kdist[-1])
    ax_band.set_ylim(energy_range)
    ax_band.grid(True, alpha=0.3, linestyle='--', axis='y')
    ax_band.tick_params(labelsize=12)
    ax_band.legend(fontsize=11, loc='upper right')

    # Right panel: PDOS
    ax_dos = fig.add_subplot(gs[1], sharey=ax_band)

    # Plot PDOS for each group
    colors = plt.cm.get_cmap('tab10')(np.linspace(0, 1, len(pdos_dict)))

    for idx, (name, dos) in enumerate(pdos_dict.items()):
        ax_dos.plot(dos, dos_energy - efermi, linewidth=2.5, label=name,
                   color=colors[idx], alpha=0.8)
        ax_dos.fill_betweenx(dos_energy - efermi, 0, dos, alpha=0.3, color=colors[idx])

    # Fermi level
    ax_dos.axhline(0, color='k', linestyle='--', linewidth=1.5, alpha=0.7)

    ax_dos.set_xlabel('DOS (arb. units)', fontsize=14, fontweight='bold')
    ax_dos.set_title('PDOS', fontsize=14, fontweight='bold')
    ax_dos.set_ylim(energy_range)
    ax_dos.set_xlim(left=0)
    ax_dos.grid(True, alpha=0.3, linestyle='--', axis='y')
    ax_dos.tick_params(labelsize=12)
    ax_dos.legend(fontsize=11, loc='upper right', framealpha=0.9)
    ax_dos.yaxis.set_ticklabels([])  # Remove y-tick labels (shared with band structure)

    fig.suptitle(title, fontsize=16, fontweight='bold', y=0.98)
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig(output, dpi=300, bbox_inches='tight')
    print(f"Combined band structure and PDOS saved to: {output}")
    plt.close()


def plot_band_with_layer_pdos(
    kpoints: np.ndarray,
    energies: np.ndarray,
    layer_projections: Dict[str, np.ndarray],
    dos_energy: np.ndarray,
    pdos_dict: Dict[str, np.ndarray],
    efermi: float = 0.0,
    high_sym_points: Optional[List[int]] = None,
    high_sym_labels: Optional[List[str]] = None,
    title: str = 'Band Structure and Layer-PDOS',
    output: str = 'band_layer_pdos.png',
    energy_range: Tuple[float, float] = (-3, 3)
):
    """
    Plot band structure with layer projections + layer-PDOS

    Args:
        kpoints: k-point coordinates
        energies: Band energies
        layer_projections: Layer projections for band structure
        dos_energy: DOS energy grid
        pdos_dict: Layer-resolved PDOS
        efermi: Fermi energy
        high_sym_points: High-symmetry points
        high_sym_labels: High-symmetry labels
        title: Plot title
        output: Output filename
        energy_range: Energy window
    """
    kdist = get_kpath_distance(kpoints)
    energies_shifted = energies - efermi

    # Auto-detect high-symmetry points if not provided
    if high_sym_points is None:
        from plot_fatband import find_high_symmetry_points
        high_sym_points, high_sym_labels = find_high_symmetry_points(kpoints)

    # Create figure
    fig = plt.figure(figsize=(14, 8))
    gs = gridspec.GridSpec(1, 2, width_ratios=[2, 1], wspace=0.05)

    # Left panel: Layer-projected band structure
    ax_band = fig.add_subplot(gs[0])

    layer_names = list(layer_projections.keys())
    colors = plt.cm.get_cmap('tab10')(np.linspace(0, 1, len(layer_names)))
    layer_colors = {name: colors[i] for i, name in enumerate(layer_names)}

    # Plot bands with layer projections
    for layer_name, color in layer_colors.items():
        projections = layer_projections[layer_name]

        for iband in range(energies.shape[1]):
            band_energies = energies_shifted[:, iband]
            band_proj = projections[:, iband]

            mask = (band_energies >= energy_range[0]) & (band_energies <= energy_range[1])
            if not np.any(mask):
                continue

            # Scatter with size based on projection
            scatter_size = 30 * band_proj[mask] / projections.max()
            ax_band.scatter(
                kdist[mask],
                band_energies[mask],
                s=scatter_size,
                color=color,
                alpha=0.6,
                edgecolors='none',
                label=layer_name if iband == 0 else ''
            )

    # Fermi level
    ax_band.axhline(0, color='k', linestyle='--', linewidth=1.5, alpha=0.7)

    # High-symmetry points
    if high_sym_points is not None and high_sym_labels is not None:
        for kpt in high_sym_points:
            ax_band.axvline(kdist[kpt], color='gray', linestyle='-', linewidth=0.8, alpha=0.5)

        ax_band.set_xticks([kdist[kpt] for kpt in high_sym_points])
        ax_band.set_xticklabels(high_sym_labels, fontsize=14)

    ax_band.set_ylabel('Energy - $E_F$ (eV)', fontsize=14, fontweight='bold')
    ax_band.set_title('Layer-Resolved Band Structure', fontsize=14, fontweight='bold')
    ax_band.set_xlim(kdist[0], kdist[-1])
    ax_band.set_ylim(energy_range)
    ax_band.grid(True, alpha=0.3, linestyle='--', axis='y')
    ax_band.tick_params(labelsize=12)

    # Legend
    handles, labels = ax_band.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ax_band.legend(by_label.values(), by_label.keys(), fontsize=11,
                  loc='upper right', framealpha=0.9)

    # Right panel: Layer-PDOS
    ax_dos = fig.add_subplot(gs[1], sharey=ax_band)

    for layer_name, color in layer_colors.items():
        if layer_name in pdos_dict:
            dos = pdos_dict[layer_name]
            ax_dos.plot(dos, dos_energy - efermi, linewidth=2.5,
                       label=layer_name, color=color, alpha=0.8)
            ax_dos.fill_betweenx(dos_energy - efermi, 0, dos,
                                alpha=0.3, color=color)

    # Fermi level
    ax_dos.axhline(0, color='k', linestyle='--', linewidth=1.5, alpha=0.7)

    ax_dos.set_xlabel('PDOS (arb. units)', fontsize=14, fontweight='bold')
    ax_dos.set_title('Layer-PDOS', fontsize=14, fontweight='bold')
    ax_dos.set_ylim(energy_range)
    ax_dos.set_xlim(left=0)
    ax_dos.grid(True, alpha=0.3, linestyle='--', axis='y')
    ax_dos.tick_params(labelsize=12)
    ax_dos.legend(fontsize=11, loc='upper right', framealpha=0.9)
    ax_dos.yaxis.set_ticklabels([])

    fig.suptitle(title, fontsize=16, fontweight='bold', y=0.98)
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig(output, dpi=300, bbox_inches='tight')
    print(f"Combined band and layer-PDOS saved to: {output}")
    plt.close()


def main():
    """Main function"""
    import argparse

    parser = argparse.ArgumentParser(description='Combined band structure and PDOS plotting')
    parser.add_argument('--procar', default='PROCAR', help='PROCAR file')
    parser.add_argument('--doscar', default='DOSCAR', help='DOSCAR file')
    parser.add_argument('--output', default='band_pdos.png', help='Output file')
    parser.add_argument('--emin', type=float, default=-3.0, help='Min energy (eV)')
    parser.add_argument('--emax', type=float, default=3.0, help='Max energy (eV)')

    args = parser.parse_args()

    # Parse PROCAR
    print(f"Parsing PROCAR: {args.procar}")
    parser_obj = pp.ProcarParser(args.procar)
    data = parser_obj.parse()

    # Read DOSCAR
    print(f"Reading DOSCAR: {args.doscar}")
    dos_energy, dos_total, efermi = read_doscar_total(args.doscar)

    print(f"Fermi energy: {efermi:.4f} eV")

    # Get energies
    energies = data.energies[:, :, 0] if data.is_spin_polarized else data.energies

    # Compute PDOS (simplified - would need atom groups)
    # For now, just plot total DOS

    pdos_dict = {'Total': dos_total}

    # Plot
    plot_band_and_pdos(
        data.kpoints, energies, dos_energy, pdos_dict, efermi,
        title='Band Structure and DOS',
        output=args.output,
        energy_range=(args.emin, args.emax)
    )

    print("\nDone!")


if __name__ == '__main__':
    # Need scipy for gaussian filtering
    try:
        import scipy.ndimage
    except ImportError:
        print("Warning: scipy not available, some features may not work")

    main()
