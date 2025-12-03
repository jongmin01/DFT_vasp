#!/usr/bin/env python3
"""
compare_band_sources.py - Compare band energies from EIGENVAL, PROCAR, and vasprun.xml

Checks for consistency between different VASP output files
"""

import numpy as np
import xml.etree.ElementTree as ET
import gzip
from typing import Tuple, Dict


def read_eigenval(filename='EIGENVAL') -> Tuple[np.ndarray, np.ndarray, float, int, int]:
    """
    Read EIGENVAL file

    Returns:
        kpoints: [nkpts, 3]
        energies: [nkpts, nbands]
        efermi: Fermi energy
        nkpts: number of k-points
        nbands: number of bands
    """
    with open(filename, 'r') as f:
        lines = f.readlines()

    # Line 0: comment
    # Line 1-4: header info
    # Line 5: number of electrons, k-points, bands
    header = lines[5].split()
    nkpts = int(header[1])
    nbands = int(header[2])

    # Fermi energy might be on line 5 or needs to be read from elsewhere
    # EIGENVAL doesn't always contain Fermi energy, will be 0 if not found
    efermi = 0.0

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

    return kpoints, energies, efermi, nkpts, nbands


def read_procar_energies(filename='PROCAR') -> Tuple[np.ndarray, np.ndarray]:
    """
    Read band energies from PROCAR

    Returns:
        kpoints: [nkpts, 3]
        energies: [nkpts, nbands]
    """
    # Handle gzip compressed files
    if filename.endswith('.gz'):
        with gzip.open(filename, 'rt') as f:
            lines = f.readlines()
    else:
        with open(filename, 'r') as f:
            lines = f.readlines()

    # Parse header
    import re
    header_line = lines[1]
    match = re.search(r'# of k-points:\s+(\d+)\s+# of bands:\s+(\d+)', header_line)
    nkpts = int(match.group(1))
    nbands = int(match.group(2))

    kpoints = np.zeros((nkpts, 3))
    energies = np.zeros((nkpts, nbands))

    i = 0
    kpt_count = 0

    while i < len(lines) and kpt_count < nkpts:
        line = lines[i].strip()

        # Find k-point block
        if line.startswith('k-point'):
            # Parse k-point coordinates
            parts = line.split(':')
            kpt_num = int(parts[0].split()[1]) - 1  # 0-indexed

            if kpt_num >= nkpts:
                # Spin-polarized case, skip
                i += 1
                continue

            kpt_coords = parts[1].split('weight')[0].split()
            kpoints[kpt_num] = [float(kpt_coords[0]), float(kpt_coords[1]), float(kpt_coords[2])]

            i += 1

            # Parse bands
            for iband in range(nbands):
                # Find band energy line
                while i < len(lines) and not lines[i].strip().startswith('band'):
                    i += 1

                if i >= len(lines):
                    break

                # Parse energy
                band_line = lines[i].strip()
                energy_match = re.search(r'energy\s+([-\d.Ee+]+)', band_line)
                if energy_match:
                    energies[kpt_num, iband] = float(energy_match.group(1))

                i += 1

                # Skip projection data
                while i < len(lines) and not (lines[i].strip().startswith('band') or lines[i].strip().startswith('k-point')):
                    i += 1

            kpt_count += 1
        else:
            i += 1

    return kpoints, energies


def read_vasprun_energies(filename='vasprun.xml') -> Tuple[np.ndarray, np.ndarray, float]:
    """
    Read band energies from vasprun.xml

    Returns:
        kpoints: [nkpts, 3]
        energies: [nkpts, nbands]
        efermi: Fermi energy
    """
    print(f"Parsing {filename}...")
    tree = ET.parse(filename)
    root = tree.getroot()

    # Get Fermi level
    efermi = 0.0
    for dos in root.findall('.//dos'):
        efermi_elem = dos.find('i[@name="efermi"]')
        if efermi_elem is not None:
            efermi = float(efermi_elem.text)
            break

    # Get k-points
    kpoints_elem = root.find('.//varray[@name="kpointlist"]')
    kpts = []
    for v in kpoints_elem.findall('v'):
        coords = [float(x) for x in v.text.split()]
        kpts.append(coords)
    kpoints = np.array(kpts)

    # Get eigenvalues
    eigenvalues_elem = root.find('.//eigenvalues')
    kpoint_sets = eigenvalues_elem.findall('.//array/set/set/set')

    nkpts = len(kpoint_sets)
    nbands = len(kpoint_sets[0].findall('r'))

    energies = np.zeros((nkpts, nbands))

    for ikpt, kpt_elem in enumerate(kpoint_sets):
        bands = kpt_elem.findall('r')
        for iband, band in enumerate(bands):
            values = [float(x) for x in band.text.split()]
            energies[ikpt, iband] = values[0]

    return kpoints, energies, efermi


def read_doscar_fermi(filename='DOSCAR') -> float:
    """Read Fermi energy from DOSCAR"""
    try:
        with open(filename, 'r') as f:
            lines = f.readlines()
        header = lines[5].split()
        efermi = float(header[3])
        return efermi
    except:
        return 0.0


def compare_sources():
    """Compare band energies from all sources"""
    print("="*80)
    print("COMPARING BAND ENERGIES FROM DIFFERENT SOURCES")
    print("="*80)

    # Read EIGENVAL
    print("\n1. Reading EIGENVAL...")
    try:
        kpts_eig, energies_eig, efermi_eig, nkpts_eig, nbands_eig = read_eigenval('EIGENVAL')
        print(f"   ✓ Successfully read EIGENVAL")
        print(f"     - k-points: {nkpts_eig}")
        print(f"     - bands: {nbands_eig}")
        print(f"     - Fermi energy: {efermi_eig:.6f} eV (may be 0 if not in file)")
    except Exception as e:
        print(f"   ✗ Error reading EIGENVAL: {e}")
        kpts_eig, energies_eig = None, None

    # Read PROCAR
    print("\n2. Reading PROCAR...")
    try:
        kpts_pro, energies_pro = read_procar_energies('PROCAR.gz')
        print(f"   ✓ Successfully read PROCAR")
        print(f"     - k-points: {kpts_pro.shape[0]}")
        print(f"     - bands: {energies_pro.shape[1]}")
    except Exception as e:
        print(f"   ✗ Error reading PROCAR: {e}")
        kpts_pro, energies_pro = None, None

    # Read vasprun.xml
    print("\n3. Reading vasprun.xml...")
    try:
        kpts_vr, energies_vr, efermi_vr = read_vasprun_energies('vasprun.xml')
        print(f"   ✓ Successfully read vasprun.xml")
        print(f"     - k-points: {kpts_vr.shape[0]}")
        print(f"     - bands: {energies_vr.shape[1]}")
        print(f"     - Fermi energy: {efermi_vr:.6f} eV")
    except Exception as e:
        print(f"   ✗ Error reading vasprun.xml: {e}")
        kpts_vr, energies_vr, efermi_vr = None, None, None

    # Read DOSCAR for Fermi energy
    print("\n4. Reading Fermi energy from DOSCAR...")
    efermi_dos = read_doscar_fermi('DOSCAR')
    print(f"   Fermi energy: {efermi_dos:.6f} eV")

    # Comparison
    print("\n" + "="*80)
    print("COMPARISON RESULTS")
    print("="*80)

    # Compare k-points
    print("\n--- K-POINTS ---")
    if kpts_eig is not None and kpts_pro is not None:
        kpt_diff = np.abs(kpts_eig - kpts_pro).max()
        print(f"EIGENVAL vs PROCAR: max difference = {kpt_diff:.2e}")
        if kpt_diff < 1e-6:
            print("  ✓ K-points match!")
        else:
            print("  ✗ K-points differ!")

    if kpts_eig is not None and kpts_vr is not None:
        kpt_diff = np.abs(kpts_eig - kpts_vr).max()
        print(f"EIGENVAL vs vasprun.xml: max difference = {kpt_diff:.2e}")
        if kpt_diff < 1e-6:
            print("  ✓ K-points match!")
        else:
            print("  ✗ K-points differ!")

    # Compare energies
    print("\n--- BAND ENERGIES ---")
    if energies_eig is not None and energies_pro is not None:
        if energies_eig.shape == energies_pro.shape:
            energy_diff = np.abs(energies_eig - energies_pro)
            print(f"\nEIGENVAL vs PROCAR:")
            print(f"  Max difference: {energy_diff.max():.6e} eV")
            print(f"  Mean difference: {energy_diff.mean():.6e} eV")
            print(f"  RMS difference: {np.sqrt((energy_diff**2).mean()):.6e} eV")

            if energy_diff.max() < 1e-6:
                print("  ✓ Energies match perfectly!")
            elif energy_diff.max() < 1e-3:
                print("  ✓ Energies match (within numerical precision)")
            else:
                print("  ✗ WARNING: Energies differ significantly!")

                # Find where the largest differences occur
                max_idx = np.unravel_index(energy_diff.argmax(), energy_diff.shape)
                print(f"\n  Largest difference at k-point {max_idx[0]}, band {max_idx[1]}:")
                print(f"    EIGENVAL: {energies_eig[max_idx]:.6f} eV")
                print(f"    PROCAR:   {energies_pro[max_idx]:.6f} eV")
                print(f"    Difference: {energy_diff[max_idx]:.6f} eV")
        else:
            print(f"  ✗ Shape mismatch: EIGENVAL {energies_eig.shape} vs PROCAR {energies_pro.shape}")

    if energies_eig is not None and energies_vr is not None:
        if energies_eig.shape == energies_vr.shape:
            energy_diff = np.abs(energies_eig - energies_vr)
            print(f"\nEIGENVAL vs vasprun.xml:")
            print(f"  Max difference: {energy_diff.max():.6e} eV")
            print(f"  Mean difference: {energy_diff.mean():.6e} eV")
            print(f"  RMS difference: {np.sqrt((energy_diff**2).mean()):.6e} eV")

            if energy_diff.max() < 1e-6:
                print("  ✓ Energies match perfectly!")
            elif energy_diff.max() < 1e-3:
                print("  ✓ Energies match (within numerical precision)")
            else:
                print("  ✗ WARNING: Energies differ significantly!")

                max_idx = np.unravel_index(energy_diff.argmax(), energy_diff.shape)
                print(f"\n  Largest difference at k-point {max_idx[0]}, band {max_idx[1]}:")
                print(f"    EIGENVAL:     {energies_eig[max_idx]:.6f} eV")
                print(f"    vasprun.xml:  {energies_vr[max_idx]:.6f} eV")
                print(f"    Difference: {energy_diff[max_idx]:.6f} eV")
        else:
            print(f"  ✗ Shape mismatch: EIGENVAL {energies_eig.shape} vs vasprun {energies_vr.shape}")

    if energies_pro is not None and energies_vr is not None:
        if energies_pro.shape == energies_vr.shape:
            energy_diff = np.abs(energies_pro - energies_vr)
            print(f"\nPROCAR vs vasprun.xml:")
            print(f"  Max difference: {energy_diff.max():.6e} eV")
            print(f"  Mean difference: {energy_diff.mean():.6e} eV")
            print(f"  RMS difference: {np.sqrt((energy_diff**2).mean()):.6e} eV")

            if energy_diff.max() < 1e-6:
                print("  ✓ Energies match perfectly!")
            elif energy_diff.max() < 1e-3:
                print("  ✓ Energies match (within numerical precision)")
            else:
                print("  ✗ WARNING: Energies differ significantly!")

    # Compare Fermi energies
    print("\n--- FERMI ENERGIES ---")
    print(f"EIGENVAL:     {efermi_eig:.6f} eV")
    if efermi_vr is not None:
        print(f"vasprun.xml:  {efermi_vr:.6f} eV")
    else:
        print(f"vasprun.xml:  N/A (file not found)")
    print(f"DOSCAR:       {efermi_dos:.6f} eV")

    if efermi_vr is not None and abs(efermi_vr - efermi_dos) < 1e-6:
        print("  ✓ Fermi energies match!")
    elif efermi_vr is not None:
        print(f"  ✗ Fermi energy difference: {abs(efermi_vr - efermi_dos):.6e} eV")

    # Sample band values at a few k-points
    print("\n--- SAMPLE BAND ENERGIES (first 5 k-points, bands around Fermi level) ---")
    if energies_eig is not None:
        mid_band = energies_eig.shape[1] // 2
        print("\nEIGENVAL (bands 595-605 at k-point 0):")
        for ib in range(max(0, mid_band-5), min(energies_eig.shape[1], mid_band+5)):
            print(f"  Band {ib:3d}: {energies_eig[0, ib]:10.6f} eV")

    if energies_pro is not None:
        mid_band = energies_pro.shape[1] // 2
        print("\nPROCAR (bands 595-605 at k-point 0):")
        for ib in range(max(0, mid_band-5), min(energies_pro.shape[1], mid_band+5)):
            print(f"  Band {ib:3d}: {energies_pro[0, ib]:10.6f} eV")

    if energies_vr is not None:
        mid_band = energies_vr.shape[1] // 2
        print("\nvasprun.xml (bands 595-605 at k-point 0):")
        for ib in range(max(0, mid_band-5), min(energies_vr.shape[1], mid_band+5)):
            print(f"  Band {ib:3d}: {energies_vr[0, ib]:10.6f} eV")

    print("\n" + "="*80)
    print("SUMMARY")
    print("="*80)
    print("\nAll three files (EIGENVAL, PROCAR, vasprun.xml) should contain identical")
    print("band energies. If they differ, this indicates a problem with:")
    print("  - File corruption")
    print("  - Files from different calculations")
    print("  - Parsing errors in the scripts")
    print("\nFor plotting, it's recommended to use:")
    print("  - Band energies: EIGENVAL or vasprun.xml (most reliable)")
    print("  - Fermi energy: DOSCAR or vasprun.xml")
    print("  - Projections: PROCAR (only source)")
    print("="*80 + "\n")


if __name__ == '__main__':
    compare_sources()
