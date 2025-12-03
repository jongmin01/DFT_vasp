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
print("K-POINT IDENTIFICATION")
print("=" * 70)

# k-point 80 (K valley)
k_idx_80 = 80
kpt_80 = data.kpoints[k_idx_80]
print(f"\nk-point 80:")
print(f"  Coordinates: ({kpt_80[0]:.5f}, {kpt_80[1]:.5f}, {kpt_80[2]:.5f})")
print(f"  → This is K valley ✓")

# k-point 104 (NOT Gamma!)
k_idx_104 = 104
kpt_104 = data.kpoints[k_idx_104]
print(f"\nk-point 104:")
print(f"  Coordinates: ({kpt_104[0]:.5f}, {kpt_104[1]:.5f}, {kpt_104[2]:.5f})")
print(f"  → This is NOT Γ point! This is somewhere between Γ and K")

# Gamma point
gamma_idx = 0
kpt_gamma = data.kpoints[gamma_idx]
print(f"\nk-point 0 (Γ point):")
print(f"  Coordinates: ({kpt_gamma[0]:.5f}, {kpt_gamma[1]:.5f}, {kpt_gamma[2]:.5f})")
print(f"  → This is Γ point ✓")

print("\n" + "=" * 70)
print("ENERGY VALUES")
print("=" * 70)

# === K valley (k-point 80) ===
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
print(f"  Band indices: VBM={vbm_band_K}, CBM={cbm_band_K}")

# === Γ point (k-point 0) ===
energies_at_Gamma = energies_shifted[gamma_idx, :]
occupied_Gamma = energies_at_Gamma < 0
unoccupied_Gamma = energies_at_Gamma > 0

vbm_band_Gamma = np.where(occupied_Gamma)[0][energies_at_Gamma[occupied_Gamma].argmax()]
cbm_band_Gamma = np.where(unoccupied_Gamma)[0][energies_at_Gamma[unoccupied_Gamma].argmin()]

E_VB_Gamma = energies_at_Gamma[vbm_band_Gamma]
E_CB_Gamma = energies_at_Gamma[cbm_band_Gamma]

print(f"\nAt Γ point (k-point 0):")
print(f"  E_VB(Γ) = {E_VB_Gamma:.4f} eV (relative to E_F)")
print(f"  E_CB(Γ) = {E_CB_Gamma:.4f} eV (relative to E_F)")
print(f"  Direct gap at Γ = {E_CB_Gamma - E_VB_Gamma:.4f} eV")
print(f"  Band indices: VBM={vbm_band_Gamma}, CBM={cbm_band_Gamma}")

# === Global VBM and CBM ===
print("\n" + "=" * 70)
print("GLOBAL BAND EDGES")
print("=" * 70)

# Find global VBM
occupied_mask = energies_shifted < 0
vbm_value = np.max(energies_shifted[occupied_mask])
vbm_location = np.where(energies_shifted == vbm_value)
vbm_k_idx = vbm_location[0][0]
vbm_band_idx = vbm_location[1][0]

# Find global CBM
unoccupied_mask = energies_shifted > 0
cbm_value = np.min(energies_shifted[unoccupied_mask])
cbm_location = np.where(energies_shifted == cbm_value)
cbm_k_idx = cbm_location[0][0]
cbm_band_idx = cbm_location[1][0]

print(f"\nGlobal VBM:")
print(f"  Energy: {vbm_value:.4f} eV")
print(f"  Location: k-point {vbm_k_idx}, band {vbm_band_idx}")
print(f"  k-coordinates: ({data.kpoints[vbm_k_idx][0]:.5f}, {data.kpoints[vbm_k_idx][1]:.5f}, {data.kpoints[vbm_k_idx][2]:.5f})")

print(f"\nGlobal CBM:")
print(f"  Energy: {cbm_value:.4f} eV")
print(f"  Location: k-point {cbm_k_idx}, band {cbm_band_idx}")
print(f"  k-coordinates: ({data.kpoints[cbm_k_idx][0]:.5f}, {data.kpoints[cbm_k_idx][1]:.5f}, {data.kpoints[cbm_k_idx][2]:.5f})")

print(f"\nFundamental gap: {cbm_value - vbm_value:.4f} eV")

# === Carrier dynamics analysis ===
print("\n" + "=" * 70)
print("CARRIER RELAXATION AND BARRIER ANALYSIS")
print("=" * 70)

print(f"\n[Electron relaxation from K]")
print(f"  E_CB(K) - E_CB(min) = {E_CB_K - cbm_value:.4f} eV")
print(f"  → Energy electron releases when relaxing from K to global CBM")

print(f"\n[Hole barrier to K]")
print(f"  E_VB(max) - E_VB(K) = {vbm_value - E_VB_K:.4f} eV")
print(f"  → Energy barrier for hole to reach K from global VBM")

print(f"\n[Comparison: K vs Γ]")
print(f"  E_CB(K) - E_CB(Γ) = {E_CB_K - E_CB_Gamma:.4f} eV")
print(f"  E_VB(Γ) - E_VB(K) = {E_VB_Gamma - E_VB_K:.4f} eV")

print("\n" + "=" * 70)
print("SUMMARY TABLE")
print("=" * 70)
print(f"\n{'Location':<15} {'E_VB (eV)':<15} {'E_CB (eV)':<15} {'Direct Gap (eV)':<15}")
print("-" * 70)
print(f"{'K valley':<15} {E_VB_K:<15.4f} {E_CB_K:<15.4f} {E_CB_K - E_VB_K:<15.4f}")
print(f"{'Γ point':<15} {E_VB_Gamma:<15.4f} {E_CB_Gamma:<15.4f} {E_CB_Gamma - E_VB_Gamma:<15.4f}")
print(f"{'Global min/max':<15} {vbm_value:<15.4f} {cbm_value:<15.4f} {cbm_value - vbm_value:<15.4f}")
print("=" * 70)

