import procar_parser as pp
import numpy as np

# Parse PROCAR
parser = pp.ProcarParser('PROCAR.gz')
data = parser.parse()

# Read Fermi energy
def read_fermi_energy(doscar_file='DOSCAR'):
    try:
        with open(doscar_file, 'r') as f:
            lines = f.readlines()
        header = lines[5].split()
        efermi = float(header[3])
        return efermi
    except:
        return 0.0

efermi = read_fermi_energy('DOSCAR')
print(f"Fermi energy: {efermi:.4f} eV\n")

# Get energies
energies = data.energies[:, :, 0] if data.is_spin_polarized else data.energies
energies_shifted = energies - efermi

print("=" * 70)
print("CHECKING SPECIFIC K-POINTS")
print("=" * 70)

# Check k-point 80 (should be K valley)
k_idx_80 = 80
kpt_80 = data.kpoints[k_idx_80]
print(f"\nk-point 80:")
print(f"  Coordinates: ({kpt_80[0]:.5f}, {kpt_80[1]:.5f}, {kpt_80[2]:.5f})")
print(f"  Distance from K (1/3, 1/3, 0): {np.linalg.norm(kpt_80 - np.array([1/3, 1/3, 0])):.6f}")
print(f"  → This is K valley: {np.linalg.norm(kpt_80 - np.array([1/3, 1/3, 0])) < 0.01}")

# Check k-point 104
k_idx_104 = 104
kpt_104 = data.kpoints[k_idx_104]
print(f"\nk-point 104:")
print(f"  Coordinates: ({kpt_104[0]:.5f}, {kpt_104[1]:.5f}, {kpt_104[2]:.5f})")
print(f"  Distance from Γ (0, 0, 0): {np.linalg.norm(kpt_104):.6f}")
print(f"  → This is Γ point: {np.linalg.norm(kpt_104) < 0.01}")

# Find actual Gamma point
gamma_indices = []
for i, kpt in enumerate(data.kpoints):
    if np.linalg.norm(kpt) < 0.01:
        gamma_indices.append(i)

print(f"\nActual Γ points found at k-point indices: {gamma_indices}")
for idx in gamma_indices:
    kpt = data.kpoints[idx]
    print(f"  k-point {idx}: ({kpt[0]:.5f}, {kpt[1]:.5f}, {kpt[2]:.5f})")

print("\n" + "=" * 70)
print("ENERGY VALUES AT K VALLEY (k-point 80)")
print("=" * 70)

# Find VBM and CBM at K valley (k-point 80)
energies_at_K = energies_shifted[k_idx_80, :]
occupied_K = energies_at_K < 0
unoccupied_K = energies_at_K > 0

vbm_band_K = np.where(occupied_K)[0][energies_at_K[occupied_K].argmax()]
cbm_band_K = np.where(unoccupied_K)[0][energies_at_K[unoccupied_K].argmin()]

E_VB_K = energies_at_K[vbm_band_K]
E_CB_K = energies_at_K[cbm_band_K]

print(f"\nAt K valley (k-point 80):")
print(f"  E_VB(K) = {E_VB_K:.4f} eV (relative to E_F)")
print(f"  E_CB(K) = {E_CB_K:.4f} eV (relative to E_F)")
print(f"  Direct gap at K = {E_CB_K - E_VB_K:.4f} eV")
print(f"  VBM band index: {vbm_band_K}")
print(f"  CBM band index: {cbm_band_K}")

print("\n" + "=" * 70)
print("ENERGY VALUES AT Γ POINT")
print("=" * 70)

if gamma_indices:
    # Use first Gamma point (k-point 0)
    gamma_idx = gamma_indices[0]
    
    energies_at_Gamma = energies_shifted[gamma_idx, :]
    occupied_Gamma = energies_at_Gamma < 0
    unoccupied_Gamma = energies_at_Gamma > 0
    
    vbm_band_Gamma = np.where(occupied_Gamma)[0][energies_at_Gamma[occupied_Gamma].argmax()]
    cbm_band_Gamma = np.where(unoccupied_Gamma)[0][energies_at_Gamma[unoccupied_Gamma].argmin()]
    
    E_VB_Gamma = energies_at_Gamma[vbm_band_Gamma]
    E_CB_Gamma = energies_at_Gamma[cbm_band_Gamma]
    
    print(f"\nAt Γ point (k-point {gamma_idx}):")
    print(f"  E_VB(Γ) = {E_VB_Gamma:.4f} eV (relative to E_F)")
    print(f"  E_CB(Γ) = {E_CB_Gamma:.4f} eV (relative to E_F)")
    print(f"  Direct gap at Γ = {E_CB_Gamma - E_VB_Gamma:.4f} eV")
    print(f"  VBM band index: {vbm_band_Gamma}")
    print(f"  CBM band index: {cbm_band_Gamma}")

print("\n" + "=" * 70)
print("CARRIER RELAXATION AND BARRIER ANALYSIS")
print("=" * 70)

# Global VBM and CBM
all_energies_flat = energies_shifted.flatten()
occupied_all = all_energies_flat < 0
unoccupied_all = all_energies_flat > 0

global_VBM = all_energies_flat[occupied_all].max()
global_CBM = all_energies_flat[unoccupied_all].min()

# Find where global VBM and CBM occur
vbm_k_idx, vbm_band_idx = np.unravel_index(energies_shifted.argmax(where=energies_shifted<0, initial=-np.inf), energies_shifted.shape)
cbm_k_idx, cbm_band_idx = np.unravel_index(energies_shifted.argmin(where=energies_shifted>0, initial=np.inf), energies_shifted.shape)

print(f"\nGlobal band edges:")
print(f"  Global VBM = {global_VBM:.4f} eV at k-point {vbm_k_idx}")
print(f"  Global CBM = {global_CBM:.4f} eV at k-point {cbm_k_idx}")
print(f"  Fundamental gap = {global_CBM - global_VBM:.4f} eV")

print(f"\nElectron relaxation energy:")
print(f"  E_CB(K) - E_CB(min) = {E_CB_K - global_CBM:.4f} eV")
print(f"  (Energy electron loses relaxing from K to CBM)")

print(f"\nHole barrier:")
print(f"  E_VB(max) - E_VB(K) = {global_VBM - E_VB_K:.4f} eV")
print(f"  (Energy barrier for hole to reach K from VBM)")

if gamma_indices:
    print(f"\nComparison with Γ point:")
    print(f"  E_CB(K) - E_CB(Γ) = {E_CB_K - E_CB_Gamma:.4f} eV")
    print(f"  E_VB(Γ) - E_VB(K) = {E_VB_Gamma - E_VB_K:.4f} eV")

print("\n" + "=" * 70)

