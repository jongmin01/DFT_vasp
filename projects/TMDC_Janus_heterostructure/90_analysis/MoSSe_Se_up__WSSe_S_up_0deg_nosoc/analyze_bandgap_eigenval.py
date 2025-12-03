#!/usr/bin/env python3
"""
analyze_bandgap_eigenval.py - Extract Fermi level, VBM, CBM, and bandgap from EIGENVAL

Uses EIGENVAL (most reliable for band structure calculations) instead of vasprun.xml
"""

import numpy as np
from typing import Tuple, List, Dict


def read_eigenval(filename='EIGENVAL') -> Tuple[np.ndarray, np.ndarray, int, int]:
    """
    Read EIGENVAL file

    Returns:
        kpoints: [nkpts, 3]
        energies: [nkpts, nbands]
        nkpts: number of k-points
        nbands: number of bands
    """
    with open(filename, 'r') as f:
        lines = f.readlines()

    # Line 5: number of electrons, k-points, bands
    header = lines[5].split()
    nkpts = int(header[1])
    nbands = int(header[2])

    # Parse k-points and energies
    kpoints = np.zeros((nkpts, 3))
    energies = np.zeros((nkpts, nbands))

    line_idx = 6  # Start of data

    for ikpt in range(nkpts):
        # k-point line (blank line then k-point coords and weight)
        line_idx += 1  # skip blank line
        kpt_line = lines[line_idx].split()
        kpoints[ikpt] = [float(kpt_line[0]), float(kpt_line[1]), float(kpt_line[2])]
        line_idx += 1

        # Band energies
        for iband in range(nbands):
            band_line = lines[line_idx].split()
            # Format: band_index  energy  (occupancy for spin-polarized)
            energies[ikpt, iband] = float(band_line[1])
            line_idx += 1

    return kpoints, energies, nkpts, nbands


def read_doscar_fermi(filename='DOSCAR') -> float:
    """Read Fermi energy from DOSCAR"""
    try:
        with open(filename, 'r') as f:
            lines = f.readlines()
        header = lines[5].split()
        efermi = float(header[3])
        return efermi
    except:
        print("Warning: Could not read Fermi energy from DOSCAR")
        return 0.0


def find_vbm_cbm(energies: np.ndarray, efermi: float) -> Tuple[float, float, int, int, int, int]:
    """
    Find VBM and CBM based on Fermi level

    Returns:
        vbm, cbm, vbm_kpt_idx, vbm_band_idx, cbm_kpt_idx, cbm_band_idx
    """
    energies_shifted = energies - efermi

    # Find VBM: maximum energy among occupied states (E < E_F)
    occupied = energies_shifted < 0
    vbm = -np.inf
    vbm_kpt_idx = 0
    vbm_band_idx = 0

    for ikpt in range(energies.shape[0]):
        for iband in range(energies.shape[1]):
            if occupied[ikpt, iband]:
                if energies[ikpt, iband] > vbm:
                    vbm = energies[ikpt, iband]
                    vbm_kpt_idx = ikpt
                    vbm_band_idx = iband

    # Find CBM: minimum energy among unoccupied states (E > E_F)
    unoccupied = energies_shifted > 0
    cbm = np.inf
    cbm_kpt_idx = 0
    cbm_band_idx = 0

    for ikpt in range(energies.shape[0]):
        for iband in range(energies.shape[1]):
            if unoccupied[ikpt, iband]:
                if energies[ikpt, iband] < cbm:
                    cbm = energies[ikpt, iband]
                    cbm_kpt_idx = ikpt
                    cbm_band_idx = iband

    return vbm, cbm, vbm_kpt_idx, vbm_band_idx, cbm_kpt_idx, cbm_band_idx


def identify_high_symmetry_points(kpoints: np.ndarray) -> Dict[str, List[int]]:
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

    for ikpt, kpt in enumerate(kpoints):
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


def calculate_k_k_prime_gap(energies: np.ndarray, symm_points: Dict[str, List[int]],
                            efermi: float) -> Tuple[float, float]:
    """
    Calculate K-K' direct bandgap

    Args:
        energies: Band energies
        symm_points: Dictionary of high-symmetry points
        efermi: Fermi level

    Returns:
        K point gap, K' point gap (eV)
    """
    energies_shifted = energies - efermi

    k_gap = None
    kprime_gap = None

    # K point gap
    if symm_points['K']:
        ikpt = symm_points['K'][0]

        # Find VBM at this k-point
        occupied = energies_shifted[ikpt] < 0
        if np.any(occupied):
            vbm_k = energies[ikpt, occupied].max()
        else:
            vbm_k = None

        # Find CBM at this k-point
        unoccupied = energies_shifted[ikpt] > 0
        if np.any(unoccupied):
            cbm_k = energies[ikpt, unoccupied].min()
        else:
            cbm_k = None

        if vbm_k is not None and cbm_k is not None:
            k_gap = cbm_k - vbm_k

    # K' point gap
    if symm_points['K_prime']:
        ikpt = symm_points['K_prime'][0]

        occupied = energies_shifted[ikpt] < 0
        if np.any(occupied):
            vbm_kp = energies[ikpt, occupied].max()
        else:
            vbm_kp = None

        unoccupied = energies_shifted[ikpt] > 0
        if np.any(unoccupied):
            cbm_kp = energies[ikpt, unoccupied].min()
        else:
            cbm_kp = None

        if vbm_kp is not None and cbm_kp is not None:
            kprime_gap = cbm_kp - vbm_kp

    return k_gap, kprime_gap


def get_band_energies_at_point(energies: np.ndarray, kpt_idx: int, efermi: float,
                               num_bands: int = 10) -> Dict:
    """
    Get band energies around Fermi level at specific k-point

    Args:
        energies: Band energies
        kpt_idx: k-point index
        efermi: Fermi level
        num_bands: number of bands to show above and below Fermi level

    Returns:
        Dictionary with valence and conduction band energies
    """
    energies_shifted = energies[kpt_idx] - efermi

    # Find highest occupied band
    occupied = energies_shifted < 0
    if np.any(occupied):
        occupied_indices = np.where(occupied)[0]
        highest_vb = occupied_indices[-1]
        vb_start = max(0, highest_vb - num_bands + 1)
        vb_energies = energies[kpt_idx, vb_start:highest_vb+1]
    else:
        vb_energies = np.array([])

    # Find lowest unoccupied band
    unoccupied = energies_shifted > 0
    if np.any(unoccupied):
        unoccupied_indices = np.where(unoccupied)[0]
        lowest_cb = unoccupied_indices[0]
        cb_end = min(energies.shape[1], lowest_cb + num_bands)
        cb_energies = energies[kpt_idx, lowest_cb:cb_end]
    else:
        cb_energies = np.array([])

    return {
        'valence': vb_energies,
        'conduction': cb_energies
    }


def print_analysis(kpoints: np.ndarray, energies: np.ndarray, efermi: float):
    """Print comprehensive bandgap analysis"""
    print("\n" + "="*80)
    print("BANDGAP ANALYSIS (from EIGENVAL)")
    print("="*80)

    print(f"\nNumber of k-points: {len(kpoints)}")
    print(f"Number of bands: {energies.shape[1]}")

    # Fermi level
    print(f"\nFermi Level: {efermi:.6f} eV")

    # VBM and CBM
    vbm, cbm, vbm_kpt, vbm_band, cbm_kpt, cbm_band = find_vbm_cbm(energies, efermi)
    indirect_gap = cbm - vbm

    print(f"\nValence Band Maximum (VBM):")
    print(f"  Energy: {vbm:.6f} eV")
    print(f"  k-point index: {vbm_kpt}")
    print(f"  k-point coords: {kpoints[vbm_kpt]}")
    print(f"  Band index: {vbm_band}")

    print(f"\nConduction Band Minimum (CBM):")
    print(f"  Energy: {cbm:.6f} eV")
    print(f"  k-point index: {cbm_kpt}")
    print(f"  k-point coords: {kpoints[cbm_kpt]}")
    print(f"  Band index: {cbm_band}")

    print(f"\nIndirect Bandgap: {indirect_gap:.6f} eV")
    if vbm_kpt != cbm_kpt:
        print(f"  (VBM at k-point {vbm_kpt}, CBM at k-point {cbm_kpt})")
    else:
        print(f"  (Direct gap at k-point {vbm_kpt})")

    # High-symmetry points
    symm_points = identify_high_symmetry_points(kpoints)

    print(f"\n" + "-"*80)
    print("HIGH-SYMMETRY POINTS")
    print("-"*80)

    for point_name, indices in symm_points.items():
        if indices:
            print(f"\n{point_name} point(s): k-point indices {indices}")
            for idx in indices[:3]:  # Show first 3 if multiple
                energies_dict = get_band_energies_at_point(energies, idx, efermi, num_bands=5)
                print(f"  k-point {idx}: {kpoints[idx]}")
                if len(energies_dict['valence']) > 0:
                    print(f"    Top valence bands (eV): {energies_dict['valence'][-5:]}")
                if len(energies_dict['conduction']) > 0:
                    print(f"    Bottom conduction bands (eV): {energies_dict['conduction'][:5]}")

                # Calculate local gap
                if len(energies_dict['valence']) > 0 and len(energies_dict['conduction']) > 0:
                    local_gap = energies_dict['conduction'][0] - energies_dict['valence'][-1]
                    print(f"    Direct gap at this point: {local_gap:.6f} eV")

    # K-K' direct bandgap
    print(f"\n" + "-"*80)
    print("K-K' VALLEY BANDGAPS")
    print("-"*80)

    k_gap, kprime_gap = calculate_k_k_prime_gap(energies, symm_points, efermi)

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
    print(f"Fermi Level:           {efermi:.6f} eV")
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

    # Read EIGENVAL
    print("Reading EIGENVAL...")
    kpoints, energies, nkpts, nbands = read_eigenval('EIGENVAL')
    print(f"Successfully read EIGENVAL: {nkpts} k-points, {nbands} bands")

    # Read Fermi energy from DOSCAR
    print("\nReading Fermi energy from DOSCAR...")
    efermi = read_doscar_fermi('DOSCAR')
    print(f"Fermi energy: {efermi:.6f} eV")

    # Perform analysis
    print_analysis(kpoints, energies, efermi)


if __name__ == '__main__':
    main()
