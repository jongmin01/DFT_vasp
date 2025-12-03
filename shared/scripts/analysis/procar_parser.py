#!/usr/bin/env python3
"""
procar_parser.py - Comprehensive VASP PROCAR file parser

This module provides functionality to parse VASP PROCAR files and extract:
- k-point coordinates
- Band energies
- Orbital projections (s, p, d orbitals)
- Atom-resolved projections
- Spin-resolved data (if available)

Supports both plain text and gzip-compressed PROCAR files.

For TMDC Janus heterostructures (MoSSe/WSSe)
"""

import numpy as np
import re
import gzip
from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional


@dataclass
class ProcarData:
    """Container for PROCAR data"""
    nkpts: int
    nbands: int
    nions: int
    norbitals: int
    kpoints: np.ndarray  # [nkpts, 3]
    kweights: np.ndarray  # [nkpts]
    energies: np.ndarray  # [nkpts, nbands] or [nkpts, nbands, nspin]
    occupancies: np.ndarray  # [nkpts, nbands] or [nkpts, nbands, nspin]
    projections: np.ndarray  # [nkpts, nbands, nions, norbitals] or [nkpts, nbands, nions, norbitals, nspin]
    orbital_names: List[str]
    is_spin_polarized: bool
    is_soc: bool
    phase_factors: Optional[np.ndarray] = None  # For SOC calculations


class ProcarParser:
    """Parser for VASP PROCAR files"""

    def __init__(self, filename='PROCAR'):
        """
        Initialize PROCAR parser

        Args:
            filename: Path to PROCAR file
        """
        self.filename = filename
        self.data = None

    def parse(self) -> ProcarData:
        """
        Parse PROCAR file (supports .gz compressed files)

        Returns:
            ProcarData object containing all parsed data
        """
        # Check if file is gzip compressed
        if self.filename.endswith('.gz'):
            print(f"Reading compressed PROCAR file: {self.filename}")
            with gzip.open(self.filename, 'rt', encoding='utf-8') as f:
                lines = f.readlines()
        else:
            with open(self.filename, 'r') as f:
                lines = f.readlines()

        # Parse header
        header_info = self._parse_header(lines)

        # Determine calculation type
        is_spin_polarized, is_soc = self._determine_calculation_type(lines, header_info)

        # Parse k-points, bands, and projections
        kpoints, kweights, energies, occupancies, projections = self._parse_data(
            lines, header_info, is_spin_polarized, is_soc
        )

        # Create ProcarData object
        self.data = ProcarData(
            nkpts=header_info['nkpts'],
            nbands=header_info['nbands'],
            nions=header_info['nions'],
            norbitals=len(header_info['orbital_names']),
            kpoints=kpoints,
            kweights=kweights,
            energies=energies,
            occupancies=occupancies,
            projections=projections,
            orbital_names=header_info['orbital_names'],
            is_spin_polarized=is_spin_polarized,
            is_soc=is_soc
        )

        return self.data

    def _parse_header(self, lines: List[str]) -> Dict:
        """Parse PROCAR header"""
        # Line 1: "PROCAR lm decomposed" or similar
        # Line 2: "# of k-points:  X   # of bands:  Y   # of ions:  Z"

        header_line = lines[1]
        match = re.search(r'# of k-points:\s+(\d+)\s+# of bands:\s+(\d+)\s+# of ions:\s+(\d+)', header_line)

        if not match:
            raise ValueError("Could not parse PROCAR header")

        nkpts = int(match.group(1))
        nbands = int(match.group(2))
        nions = int(match.group(3))

        # Find orbital names from first ion block
        orbital_names = self._find_orbital_names(lines)

        return {
            'nkpts': nkpts,
            'nbands': nbands,
            'nions': nions,
            'orbital_names': orbital_names
        }

    def _find_orbital_names(self, lines: List[str]) -> List[str]:
        """Find orbital names from PROCAR file"""
        # Look for line with "ion s py pz px dxy dyz dz2 dxz dx2 tot"
        for i, line in enumerate(lines[:200]):  # Search in first 200 lines
            if 'ion' in line and 's' in line and 'tot' in line:
                parts = line.split()
                # Remove 'ion' and 'tot'
                orbital_names = [p for p in parts if p not in ['ion', 'tot']]
                return orbital_names

        # Default for standard VASP (l=0,1,2)
        return ['s', 'py', 'pz', 'px', 'dxy', 'dyz', 'dz2', 'dxz', 'dx2']

    def _determine_calculation_type(self, lines: List[str], header_info: Dict) -> Tuple[bool, bool]:
        """Determine if calculation is spin-polarized or SOC"""
        # Check for spin polarization
        # In spin-polarized calculations, there are separate blocks for spin-up and spin-down

        # Count number of k-point blocks
        k_blocks = sum(1 for line in lines if line.strip().startswith('k-point'))

        expected_blocks_nospin = header_info['nkpts']
        expected_blocks_spin = header_info['nkpts'] * 2

        is_spin_polarized = (k_blocks == expected_blocks_spin)

        # Check for SOC (phase factors present)
        is_soc = any('phase' in line.lower() for line in lines[:1000])

        return is_spin_polarized, is_soc

    def _parse_data(self, lines: List[str], header_info: Dict,
                    is_spin_polarized: bool, is_soc: bool) -> Tuple:
        """Parse k-points, energies, and projections"""
        nkpts = header_info['nkpts']
        nbands = header_info['nbands']
        nions = header_info['nions']
        norbitals = len(header_info['orbital_names'])

        # Initialize arrays
        if is_spin_polarized:
            energies = np.zeros((nkpts, nbands, 2))
            occupancies = np.zeros((nkpts, nbands, 2))
            projections = np.zeros((nkpts, nbands, nions, norbitals, 2))
        else:
            energies = np.zeros((nkpts, nbands))
            occupancies = np.zeros((nkpts, nbands))
            projections = np.zeros((nkpts, nbands, nions, norbitals))

        kpoints = np.zeros((nkpts, 3))
        kweights = np.zeros(nkpts)

        # Parse data
        i = 0
        while i < len(lines):
            line = lines[i].strip()

            # Find k-point block
            if line.startswith('k-point'):
                # Parse k-point info
                # Format: "k-point    1 :    0.00000000 0.00000000 0.00000000     weight = 0.00800000"
                parts = line.split(':')
                kpt_num = int(parts[0].split()[1]) - 1  # 0-indexed

                if is_spin_polarized and kpt_num >= nkpts:
                    # This is spin-down block
                    kpt_num = kpt_num - nkpts
                    spin_idx = 1
                else:
                    spin_idx = 0

                kpt_coords = parts[1].split('weight')[0].split()
                kpoints[kpt_num] = [float(kpt_coords[0]), float(kpt_coords[1]), float(kpt_coords[2])]

                weight_match = re.search(r'weight\s*=\s*([\d.Ee+-]+)', line)
                if weight_match:
                    kweights[kpt_num] = float(weight_match.group(1))

                i += 1

                # Parse bands for this k-point
                for iband in range(nbands):
                    # Find band energy line
                    while i < len(lines) and not lines[i].strip().startswith('band'):
                        i += 1

                    if i >= len(lines):
                        break

                    # Parse band energy
                    # Format: "band   1 # energy   -38.71370000 # occ.  1.00000000"
                    band_line = lines[i].strip()
                    energy_match = re.search(r'energy\s+([-\d.Ee+]+)', band_line)
                    occ_match = re.search(r'occ\.\s+([-\d.Ee+]+)', band_line)

                    if energy_match:
                        if is_spin_polarized:
                            energies[kpt_num, iband, spin_idx] = float(energy_match.group(1))
                        else:
                            energies[kpt_num, iband] = float(energy_match.group(1))

                    if occ_match:
                        if is_spin_polarized:
                            occupancies[kpt_num, iband, spin_idx] = float(occ_match.group(1))
                        else:
                            occupancies[kpt_num, iband] = float(occ_match.group(1))

                    i += 1

                    # Skip to ion data (skip header lines)
                    while i < len(lines) and 'ion' not in lines[i]:
                        i += 1
                    i += 1  # Skip the "ion s py pz..." header

                    # Parse ion projections
                    for iion in range(nions):
                        if i >= len(lines):
                            break

                        ion_line = lines[i].strip().split()

                        # Format: "1  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.000"
                        # ion_num, s, py, pz, px, dxy, dyz, dz2, dxz, dx2, tot

                        if len(ion_line) >= norbitals + 1:
                            for iorb in range(norbitals):
                                if is_spin_polarized:
                                    projections[kpt_num, iband, iion, iorb, spin_idx] = float(ion_line[iorb + 1])
                                else:
                                    projections[kpt_num, iband, iion, iorb] = float(ion_line[iorb + 1])

                        i += 1

                    # Skip tot line
                    i += 1

            i += 1

        return kpoints, kweights, energies, occupancies, projections

    def get_orbital_projection(self, orbital_name: str) -> np.ndarray:
        """
        Get projection for specific orbital

        Args:
            orbital_name: Name of orbital (e.g., 's', 'px', 'dz2')

        Returns:
            Orbital projection array [nkpts, nbands, nions] or [nkpts, nbands, nions, nspin]
        """
        if self.data is None:
            raise ValueError("No data loaded. Call parse() first.")

        if orbital_name not in self.data.orbital_names:
            raise ValueError(f"Orbital '{orbital_name}' not found. Available: {self.data.orbital_names}")

        orb_idx = self.data.orbital_names.index(orbital_name)

        if self.data.is_spin_polarized:
            return self.data.projections[:, :, :, orb_idx, :]
        else:
            return self.data.projections[:, :, :, orb_idx]

    def get_atom_projection(self, atom_indices: List[int]) -> np.ndarray:
        """
        Get projection for specific atoms (summed over orbitals)

        Args:
            atom_indices: List of atom indices (0-indexed)

        Returns:
            Atom projection array [nkpts, nbands] or [nkpts, nbands, nspin]
        """
        if self.data is None:
            raise ValueError("No data loaded. Call parse() first.")

        # Sum over specified atoms and all orbitals
        if self.data.is_spin_polarized:
            proj = self.data.projections[:, :, atom_indices, :, :].sum(axis=(2, 3))
        else:
            proj = self.data.projections[:, :, atom_indices, :].sum(axis=(2, 3))

        return proj

    def get_layer_projection(self, layer_atoms: Dict[str, List[int]]) -> Dict[str, np.ndarray]:
        """
        Get projections for different layers

        Args:
            layer_atoms: Dictionary mapping layer names to atom indices
                        e.g., {'MoSSe': [0,1,2], 'WSSe': [3,4,5]}

        Returns:
            Dictionary of layer projections
        """
        layer_proj = {}

        for layer_name, atom_indices in layer_atoms.items():
            layer_proj[layer_name] = self.get_atom_projection(atom_indices)

        return layer_proj


def read_poscar_for_atoms(poscar_file: str) -> Dict[str, List[int]]:
    """
    Read POSCAR to get atom types and indices

    Args:
        poscar_file: Path to POSCAR file

    Returns:
        Dictionary with atom information
    """
    with open(poscar_file, 'r') as f:
        lines = f.readlines()

    # Line 5: atom types
    # Line 6: number of each atom type

    atom_types = lines[5].split()
    atom_counts = [int(x) for x in lines[6].split()]

    # Create atom index mapping
    atom_info = {}
    current_idx = 0

    for atom_type, count in zip(atom_types, atom_counts):
        indices = list(range(current_idx, current_idx + count))
        atom_info[atom_type] = indices
        current_idx += count

    return atom_info


if __name__ == '__main__':
    import sys

    if len(sys.argv) > 1:
        procar_file = sys.argv[1]
    else:
        procar_file = 'PROCAR'

    print(f"Parsing PROCAR file: {procar_file}")

    parser = ProcarParser(procar_file)
    data = parser.parse()

    print(f"\nPROCAR Summary:")
    print(f"  Number of k-points: {data.nkpts}")
    print(f"  Number of bands: {data.nbands}")
    print(f"  Number of ions: {data.nions}")
    print(f"  Number of orbitals: {data.norbitals}")
    print(f"  Orbitals: {data.orbital_names}")
    print(f"  Spin-polarized: {data.is_spin_polarized}")
    print(f"  SOC: {data.is_soc}")
    print(f"\nEnergy range: {data.energies.min():.3f} to {data.energies.max():.3f} eV")
