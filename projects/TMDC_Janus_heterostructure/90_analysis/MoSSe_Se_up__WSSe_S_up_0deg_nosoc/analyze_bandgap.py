#!/usr/bin/env python3
"""
analyze_bandgap.py - Extract Fermi level, VBM, CBM, and bandgap information from vasprun.xml

For TMDC Janus heterostructures - analyzes indirect and K-K' direct bandgaps
"""

import xml.etree.ElementTree as ET
import numpy as np
from typing import Dict, Tuple, List, Optional


class BandgapAnalyzer:
    """Analyze bandgap from vasprun.xml"""

    def __init__(self, vasprun_file='vasprun.xml'):
        """
        Initialize analyzer

        Args:
            vasprun_file: Path to vasprun.xml
        """
        self.vasprun_file = vasprun_file
        self.efermi = None
        self.kpoints = None
        self.kpoint_coords = None
        self.eigenvalues = None
        self.occupancies = None
        self.nbands = None
        self.nkpts = None

    def parse_vasprun(self):
        """Parse vasprun.xml file"""
        print(f"Parsing {self.vasprun_file}...")
        print("This may take a while for large files...")

        tree = ET.parse(self.vasprun_file)
        root = tree.getroot()

        # Get Fermi level
        self._extract_fermi_level(root)

        # Get k-points
        self._extract_kpoints(root)

        # Get eigenvalues
        self._extract_eigenvalues(root)

        print(f"Parsed successfully!")
        print(f"  Fermi level: {self.efermi:.4f} eV")
        print(f"  Number of k-points: {self.nkpts}")
        print(f"  Number of bands: {self.nbands}")

    def _extract_fermi_level(self, root):
        """Extract Fermi level from vasprun.xml"""
        # Find efermi in DOS section
        for dos in root.findall('.//dos'):
            efermi_elem = dos.find('i[@name="efermi"]')
            if efermi_elem is not None:
                self.efermi = float(efermi_elem.text)
                break

    def _extract_kpoints(self, root):
        """Extract k-points from vasprun.xml"""
        # Find k-points in kpoints section
        kpoints_elem = root.find('.//varray[@name="kpointlist"]')
        if kpoints_elem is not None:
            kpts = []
            for v in kpoints_elem.findall('v'):
                coords = [float(x) for x in v.text.split()]
                kpts.append(coords)
            self.kpoint_coords = np.array(kpts)
            self.nkpts = len(kpts)

    def _extract_eigenvalues(self, root):
        """Extract eigenvalues from vasprun.xml"""
        # Find eigenvalues in calculation section
        eigenvalues_elem = root.find('.//eigenvalues')
        if eigenvalues_elem is not None:
            array_elem = eigenvalues_elem.find('.//array/set/set/set')
            if array_elem is not None:
                # Get dimensions from first k-point
                kpoint_sets = eigenvalues_elem.findall('.//array/set/set/set')
                self.nkpts = len(kpoint_sets)

                # Get number of bands from first k-point
                first_kpt = kpoint_sets[0]
                bands = first_kpt.findall('r')
                self.nbands = len(bands)

                # Initialize arrays
                self.eigenvalues = np.zeros((self.nkpts, self.nbands))
                self.occupancies = np.zeros((self.nkpts, self.nbands))

                # Parse all eigenvalues
                for ikpt, kpt_elem in enumerate(kpoint_sets):
                    bands = kpt_elem.findall('r')
                    for iband, band in enumerate(bands):
                        values = [float(x) for x in band.text.split()]
                        self.eigenvalues[ikpt, iband] = values[0]  # Energy
                        self.occupancies[ikpt, iband] = values[1]  # Occupancy

    def find_vbm_cbm(self) -> Tuple[float, float, int, int, int, int]:
        """
        Find VBM (Valence Band Maximum) and CBM (Conduction Band Minimum)

        Returns:
            vbm, cbm, vbm_kpt_idx, vbm_band_idx, cbm_kpt_idx, cbm_band_idx
        """
        # Find highest occupied band
        occupied = self.occupancies > 0.5

        # VBM: maximum energy among occupied states
        vbm = -np.inf
        vbm_kpt_idx = 0
        vbm_band_idx = 0

        for ikpt in range(self.nkpts):
            for iband in range(self.nbands):
                if occupied[ikpt, iband]:
                    if self.eigenvalues[ikpt, iband] > vbm:
                        vbm = self.eigenvalues[ikpt, iband]
                        vbm_kpt_idx = ikpt
                        vbm_band_idx = iband

        # CBM: minimum energy among unoccupied states
        cbm = np.inf
        cbm_kpt_idx = 0
        cbm_band_idx = 0

        for ikpt in range(self.nkpts):
            for iband in range(self.nbands):
                if not occupied[ikpt, iband]:
                    if self.eigenvalues[ikpt, iband] < cbm:
                        cbm = self.eigenvalues[ikpt, iband]
                        cbm_kpt_idx = ikpt
                        cbm_band_idx = iband

        return vbm, cbm, vbm_kpt_idx, vbm_band_idx, cbm_kpt_idx, cbm_band_idx

    def identify_high_symmetry_points(self) -> Dict[str, List[int]]:
        """
        Identify high-symmetry points (Gamma, K, M) from k-point coordinates

        Returns:
            Dictionary mapping point names to k-point indices
        """
        symm_points = {
            'Gamma': [],
            'K': [],
            'K_prime': [],
            'M': []
        }

        tolerance = 1e-3

        for ikpt, kpt in enumerate(self.kpoint_coords):
            # Gamma point: (0, 0, 0)
            if np.allclose(kpt, [0, 0, 0], atol=tolerance):
                symm_points['Gamma'].append(ikpt)

            # K point: (1/3, 1/3, 0) in fractional coordinates
            elif np.allclose(kpt, [1/3, 1/3, 0], atol=tolerance):
                symm_points['K'].append(ikpt)

            # K' point: (2/3, 2/3, 0) or (-1/3, -1/3, 0)
            elif (np.allclose(kpt, [2/3, 2/3, 0], atol=tolerance) or
                  np.allclose(kpt, [-1/3, -1/3, 0], atol=tolerance)):
                symm_points['K_prime'].append(ikpt)

            # M point: (0.5, 0, 0) or (0, 0.5, 0) or (0.5, 0.5, 0)
            elif (np.allclose(kpt, [0.5, 0, 0], atol=tolerance) or
                  np.allclose(kpt, [0, 0.5, 0], atol=tolerance) or
                  np.allclose(kpt, [0.5, 0.5, 0], atol=tolerance)):
                symm_points['M'].append(ikpt)

        return symm_points

    def calculate_k_k_prime_gap(self, symm_points: Dict[str, List[int]]) -> Tuple[float, float]:
        """
        Calculate K-K' direct bandgap

        Args:
            symm_points: Dictionary of high-symmetry points

        Returns:
            K point gap, K' point gap (eV)
        """
        # Find band edges at K points
        k_gap = None
        kprime_gap = None

        occupied = self.occupancies > 0.5

        # K point gap
        if symm_points['K']:
            ikpt = symm_points['K'][0]

            # Find VBM at this k-point
            vbm_k = -np.inf
            for iband in range(self.nbands):
                if occupied[ikpt, iband]:
                    vbm_k = max(vbm_k, self.eigenvalues[ikpt, iband])

            # Find CBM at this k-point
            cbm_k = np.inf
            for iband in range(self.nbands):
                if not occupied[ikpt, iband]:
                    cbm_k = min(cbm_k, self.eigenvalues[ikpt, iband])

            k_gap = cbm_k - vbm_k

        # K' point gap
        if symm_points['K_prime']:
            ikpt = symm_points['K_prime'][0]

            # Find VBM at this k-point
            vbm_kp = -np.inf
            for iband in range(self.nbands):
                if occupied[ikpt, iband]:
                    vbm_kp = max(vbm_kp, self.eigenvalues[ikpt, iband])

            # Find CBM at this k-point
            cbm_kp = np.inf
            for iband in range(self.nbands):
                if not occupied[ikpt, iband]:
                    cbm_kp = min(cbm_kp, self.eigenvalues[ikpt, iband])

            kprime_gap = cbm_kp - vbm_kp

        return k_gap, kprime_gap

    def get_band_energies_at_point(self, kpt_idx: int, num_bands: int = 10) -> Dict[str, np.ndarray]:
        """
        Get band energies around Fermi level at specific k-point

        Args:
            kpt_idx: k-point index
            num_bands: number of bands to show above and below Fermi level

        Returns:
            Dictionary with valence and conduction band energies
        """
        occupied = self.occupancies[kpt_idx] > 0.5

        # Find highest occupied band
        vb_indices = np.where(occupied)[0]
        if len(vb_indices) > 0:
            highest_vb = vb_indices[-1]
            vb_start = max(0, highest_vb - num_bands + 1)
            vb_energies = self.eigenvalues[kpt_idx, vb_start:highest_vb+1]
        else:
            vb_energies = np.array([])

        # Find lowest unoccupied band
        cb_indices = np.where(~occupied)[0]
        if len(cb_indices) > 0:
            lowest_cb = cb_indices[0]
            cb_end = min(self.nbands, lowest_cb + num_bands)
            cb_energies = self.eigenvalues[kpt_idx, lowest_cb:cb_end]
        else:
            cb_energies = np.array([])

        return {
            'valence': vb_energies,
            'conduction': cb_energies,
            'kpoint': self.kpoint_coords[kpt_idx]
        }

    def print_analysis(self):
        """Print comprehensive bandgap analysis"""
        print("\n" + "="*80)
        print("BANDGAP ANALYSIS")
        print("="*80)

        # Fermi level
        print(f"\nFermi Level: {self.efermi:.6f} eV")

        # VBM and CBM
        vbm, cbm, vbm_kpt, vbm_band, cbm_kpt, cbm_band = self.find_vbm_cbm()
        indirect_gap = cbm - vbm

        print(f"\nValence Band Maximum (VBM):")
        print(f"  Energy: {vbm:.6f} eV")
        print(f"  k-point index: {vbm_kpt}")
        print(f"  k-point coords: {self.kpoint_coords[vbm_kpt]}")
        print(f"  Band index: {vbm_band}")

        print(f"\nConduction Band Minimum (CBM):")
        print(f"  Energy: {cbm:.6f} eV")
        print(f"  k-point index: {cbm_kpt}")
        print(f"  k-point coords: {self.kpoint_coords[cbm_kpt]}")
        print(f"  Band index: {cbm_band}")

        print(f"\nIndirect Bandgap: {indirect_gap:.6f} eV")
        if vbm_kpt != cbm_kpt:
            print(f"  (VBM at k-point {vbm_kpt}, CBM at k-point {cbm_kpt})")
        else:
            print(f"  (Direct gap at k-point {vbm_kpt})")

        # High-symmetry points
        symm_points = self.identify_high_symmetry_points()

        print(f"\n" + "-"*80)
        print("HIGH-SYMMETRY POINTS")
        print("-"*80)

        for point_name, indices in symm_points.items():
            if indices:
                print(f"\n{point_name} point(s): k-point indices {indices}")
                for idx in indices[:3]:  # Show first 3 if multiple
                    energies = self.get_band_energies_at_point(idx, num_bands=5)
                    print(f"  k-point {idx}: {self.kpoint_coords[idx]}")
                    print(f"    Top valence bands (eV): {energies['valence'][-5:]}")
                    print(f"    Bottom conduction bands (eV): {energies['conduction'][:5]}")

                    # Calculate local gap
                    if len(energies['valence']) > 0 and len(energies['conduction']) > 0:
                        local_gap = energies['conduction'][0] - energies['valence'][-1]
                        print(f"    Direct gap at this point: {local_gap:.6f} eV")

        # K-K' direct bandgap
        print(f"\n" + "-"*80)
        print("K-K' VALLEY BANDGAPS")
        print("-"*80)

        k_gap, kprime_gap = self.calculate_k_k_prime_gap(symm_points)

        if k_gap is not None:
            print(f"\nK valley direct bandgap: {k_gap:.6f} eV")
        else:
            print(f"\nK valley: Not found in k-point mesh")

        if kprime_gap is not None:
            print(f"K' valley direct bandgap: {kprime_gap:.6f} eV")
        else:
            print(f"K' valley: Not found in k-point mesh")

        # Summary
        print(f"\n" + "="*80)
        print("SUMMARY")
        print("="*80)
        print(f"Fermi Level:           {self.efermi:.6f} eV")
        print(f"VBM:                   {vbm:.6f} eV (k-point {vbm_kpt})")
        print(f"CBM:                   {cbm:.6f} eV (k-point {cbm_kpt})")
        print(f"Indirect Bandgap:      {indirect_gap:.6f} eV")
        if k_gap is not None:
            print(f"K valley direct gap:   {k_gap:.6f} eV")
        if kprime_gap is not None:
            print(f"K' valley direct gap:  {kprime_gap:.6f} eV")
        print("="*80 + "\n")


def main():
    """Main function"""
    import sys

    if len(sys.argv) > 1:
        vasprun_file = sys.argv[1]
    else:
        vasprun_file = 'vasprun.xml'

    analyzer = BandgapAnalyzer(vasprun_file)
    analyzer.parse_vasprun()
    analyzer.print_analysis()


if __name__ == '__main__':
    main()
