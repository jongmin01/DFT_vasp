import os
import numpy as np
import matplotlib.pyplot as plt
from pymatgen.io.vasp.outputs import Vasprun, Procar, Eigenval
from pymatgen.electronic_structure.core import Spin
from scipy.ndimage import gaussian_filter1d

strain_dirs = ["strain099", "strain100", "strain101"]
calc_types = ["soc", "nosoc"]

results = {}

def read_doscar(doscar_file):
    with open(doscar_file, 'r') as f:
        lines = f.readlines()

    # DOS starts after 6 header lines
    data = []
    for line in lines[6:]:
        parts = line.split()
        if len(parts) < 3:
            break
        E, DOS, iDOS = map(float, parts[:3])
        data.append([E, DOS])

    data = np.array(data)
    return data[:,0], data[:,1]


for sd in strain_dirs:
    results[sd] = {}
    for calc in calc_types:
        print(f"\n===== Processing {sd}/{calc} =====")

        vasprun_gap_file = os.path.join(sd, calc, "vasprun.xml")
        
        vasprun_pdos_file = os.path.join(sd, calc, "pdos_vasprun.xml")
        eigen_file   = os.path.join(sd, calc, "EIGENVAL")

        if not (os.path.exists(vasprun_gap_file) and os.path.exists(eigen_file)):
            print(f"[!] Missing vasprun.xml/EIGENVAL in {sd}/{calc}, skipping")
            continue

        vasprun = Vasprun(vasprun_gap_file, parse_potcar_file=False) 
        structure = vasprun.final_structure
        efermi = vasprun.efermi

        # read eigenvalues to compute VBM/CBM
        eigen = Eigenval(eigen_file)
        eig_raw = eigen.eigenvalues

        try:
            occ_raw = eigen.occupancies
            if isinstance(eig_raw, dict):
                spin_key = list(eig_raw.keys())[0]
                energies = eig_raw[spin_key][:,:,0] - efermi
                occs     = occ_raw[spin_key][:,:,0]
            else:
                energies = eig_raw[:,:,0] - efermi
                occs     = occ_raw[:,:,0]
        except AttributeError:
            if isinstance(eig_raw, dict):
                spin_key = list(eig_raw.keys())[0]
                energies = eig_raw[spin_key][:,:,0] - efermi
                occs     = eig_raw[spin_key][:,:,1]
            else:
                energies = eig_raw[:,:,0] - efermi
                occs     = eig_raw[:,:,1]

        nk, nb = energies.shape

        # detect band edges
        occupied = energies[occs > 1e-3]
        empty    = energies[occs < 1e-3]
        VBM = occupied.max()
        CBM = empty.min()
        gap = CBM - VBM

        print(f"[+] {sd}/{calc}: VBM={VBM:.3f}, CBM={CBM:.3f}, gap={gap:.3f} eV")

        # ============ SOC: PDOS via PROCAR ============
        if calc == "soc":
            procar_file = os.path.join(sd, calc, "PROCAR")

            # Check for PROCAR and PDOS vasprun file
            if not os.path.exists(procar_file):
                print(f"[!] No PROCAR in SOC {sd}, skipping PDOS.")
                results[sd][calc] = {
                    "VBM": VBM, "CBM": CBM, "gap": gap,
                    "bins": None, "dz2": None,
                    "DOS_bins": None, "DOS": None
                }
                continue

            procar = Procar(procar_file)

            if isinstance(procar.data, dict):
                if Spin.up in procar.data:
                    proj = procar.data[Spin.up]
                else:
                    proj = list(procar.data.values())[0]
            else:
                proj = procar.data

            nk2, nb2, nions, norb = proj.shape
            nk_use = min(nk, nk2)
            nb_use = min(nb, nb2)

            # Mo index
            mo_idx = [i for i,s in enumerate(structure.sites) if s.species_string.lower()=="mo"]
            if len(mo_idx)==0:
                print("[!] No Mo atom found??")
                mo_idx = [0]

            # PDOS bins
            emin, emax = -3, 3
            bins = np.linspace(emin, emax, 600)
            dz2 = np.zeros_like(bins)

            for kpt in range(nk_use):
                for band in range(nb_use):
                    e = energies[kpt,band]
                    if emin < e < emax:
                        idx = np.searchsorted(bins, e) - 1
                        if 0 <= idx < len(bins):
                            dz2[idx] += proj[kpt,band,mo_idx,6].sum()

            dz2 = gaussian_filter1d(dz2, sigma=2)

            # plot PDOS
            plt.figure(figsize=(4,3), dpi=350)
            plt.plot(bins, dz2, lw=1.3, label=f"{sd}-soc")
            plt.axvline(0, color="black", linestyle="--", lw=1)
            plt.axvline(VBM, color="green", linestyle=":", lw=1)
            plt.axvline(CBM, color="red", linestyle=":", lw=1)
            plt.xlabel("Energy - $E_F$ (eV)")
            plt.ylabel("PDOS (Mo d$_{z^2}$)")
            plt.tight_layout()
            fname = f"PDOS_{sd}_SOC.png"
            plt.savefig(fname, dpi=350)
            plt.close()
            print(f"[+] saved {fname}")

            # also read total DOS for SOC
            doscar_file = os.path.join(sd, calc, "DOSCAR")
            DOS_bins, DOS_val = None, None
            if os.path.exists(doscar_file):
                Eb, Db = read_doscar(doscar_file)
                Eb -= efermi
                DOS_bins, DOS_val = Eb, gaussian_filter1d(Db, 2)

            results[sd][calc] = {
                "VBM": VBM, "CBM": CBM, "gap": gap,
                "bins": bins, "dz2": dz2,
                "DOS_bins": DOS_bins, "DOS": DOS_val
            }

        # ============ noSOC: ONLY DOS ============
        else:
            doscar_file = os.path.join(sd, calc, "DOSCAR")
            if os.path.exists(doscar_file):
                Eb, Db = read_doscar(doscar_file)
                Eb -= efermi
                DOS_val = gaussian_filter1d(Db, 2)

                plt.figure(figsize=(4,3), dpi=350)
                plt.plot(Eb, DOS_val, lw=1.3, label=f"{sd}-nosoc")
                plt.axvline(0, color="black", linestyle="--", lw=1)
                plt.xlabel("Energy - $E_F$ (eV)")
                plt.ylabel("DOS (a.u.)")
                plt.tight_layout()
                fname = f"DOS_{sd}_noSOC.png"
                plt.savefig(fname, dpi=350)
                plt.close()
                print(f"[+] saved {fname}")

                results[sd][calc] = {
                    "VBM": VBM, "CBM": CBM, "gap": gap,
                    "bins": None, "dz2": None,
                    "DOS_bins": Eb, "DOS": DOS_val
                }
            else:
                print(f"[!] no DOSCAR for nosoc {sd}")
                results[sd][calc] = {
                    "VBM": VBM, "CBM": CBM, "gap": gap,
                    "bins": None, "dz2": None,
                    "DOS_bins": None, "DOS": None
                }


# ===== Compare PDOS: SOC only =====
plt.figure(figsize=(5,3.2), dpi=400)
for sd in strain_dirs:
    if "soc" in results[sd] and results[sd]["soc"]["dz2"] is not None:
        plt.plot(results[sd]["soc"]["bins"], results[sd]["soc"]["dz2"], lw=1.3, label=f"{sd}-soc")
plt.axvline(0, color="black", linestyle="--", lw=1)
plt.xlim(-3,3)
plt.xlabel("Energy - $E_F$ (eV)")
plt.ylabel("PDOS (Mo d$_{z^2}$)")
plt.legend(frameon=False)
plt.tight_layout()
plt.savefig("compare_PDOS_SOC_strains.png", dpi=400)
plt.close()
print("[+] saved compare_PDOS_SOC_strains.png")

# ===== Compare DOS: noSOC only =====
plt.figure(figsize=(5,3.2), dpi=400)
for sd in strain_dirs:
    if "nosoc" in results[sd] and results[sd]["nosoc"]["DOS"] is not None:
        plt.plot(results[sd]["nosoc"]["DOS_bins"], results[sd]["nosoc"]["DOS"], lw=1.3, label=f"{sd}-nosoc")
plt.axvline(0, color="black", linestyle="--", lw=1)
plt.xlim(-3,3)
plt.xlabel("Energy - $E_F$ (eV)")
plt.ylabel("DOS (a.u.)")
plt.legend(frameon=False)
plt.tight_layout()
plt.savefig("compare_DOS_noSOC_strains.png", dpi=400)
plt.close()
print("[+] saved compare_DOS_noSOC_strains.png")

# ===== Same strain: SOC vs noSOC DOS =====
for sd in strain_dirs:
    if "soc" in results[sd] and "nosoc" in results[sd]:
        Eb_s = results[sd]["soc"]["DOS_bins"]
        Eb_n = results[sd]["nosoc"]["DOS_bins"]
        Db_s = results[sd]["soc"]["DOS"]
        Db_n = results[sd]["nosoc"]["DOS"]
        if Eb_s is not None and Eb_n is not None:
            plt.figure(figsize=(4.5,3), dpi=400)
            plt.plot(Eb_s, Db_s, lw=1.3, label="SOC")
            plt.plot(Eb_n, Db_n, lw=1.3, linestyle="--", label="noSOC")
            plt.axvline(0, color="black", linestyle="--", lw=1)
            plt.xlim(-3,3)
            plt.xlabel("Energy - $E_F$ (eV)")
            plt.ylabel("DOS (a.u.)")
            plt.legend(frameon=False)
            plt.title(f"{sd}: SOC vs noSOC DOS")
            plt.tight_layout()
            fname = f"compare_DOS_SOC_noSOC_{sd}.png"
            plt.savefig(fname, dpi=400)
            plt.close()
            print(f"[+] saved {fname}")

# ===== VBM / CBM / bandgap vs strain =====
str_vals = []
VBM_soc, CBM_soc, gap_soc = [], [], []
VBM_nosoc, CBM_nosoc, gap_nosoc = [], [], []

for sd in strain_dirs:
    sval = int(sd[-3:]) - 100
    str_vals.append(sval)

    if "soc" in results[sd]:
        VBM_soc.append(results[sd]["soc"]["VBM"])
        CBM_soc.append(results[sd]["soc"]["CBM"])
        gap_soc.append(results[sd]["soc"]["gap"])
    else:
        VBM_soc.append(None); CBM_soc.append(None); gap_soc.append(None)

    if "nosoc" in results[sd]:
        VBM_nosoc.append(results[sd]["nosoc"]["VBM"])
        CBM_nosoc.append(results[sd]["nosoc"]["CBM"])
        gap_nosoc.append(results[sd]["nosoc"]["gap"])
    else:
        VBM_nosoc.append(None); CBM_nosoc.append(None); gap_nosoc.append(None)

plt.figure(figsize=(5.2,3.2), dpi=400)
plt.plot(str_vals, VBM_soc, marker="o", color="green", label="VBM SOC")
plt.plot(str_vals, CBM_soc, marker="o", color="red",   label="CBM SOC")
plt.plot(str_vals, VBM_nosoc, marker="s", linestyle="--", color="green", label="VBM noSOC")
plt.plot(str_vals, CBM_nosoc, marker="s", linestyle="--", color="red",   label="CBM noSOC")
# plt.axhline(0, color="black", linestyle="--", lw=1)
plt.xlabel("Strain (%)")
plt.ylabel("Energy (eV)")
plt.legend(frameon=False, ncol=2)
plt.tight_layout()
plt.savefig("VBM_CBM_SOC_noSOC_vs_strain.png", dpi=400)
plt.close()
print("[+] saved VBM_CBM_SOC_noSOC_vs_strain.png")

plt.figure(figsize=(5.2,3.2), dpi=400)
plt.plot(str_vals, gap_soc, marker="o", color="blue", label="SOC gap")
plt.plot(str_vals, gap_nosoc, marker="s", linestyle="--", color="blue", label="noSOC gap")
plt.xlabel("Strain (%)")
plt.ylabel("Gap (eV)")
plt.legend(frameon=False)
plt.tight_layout()
plt.savefig("bandgap_SOC_noSOC_vs_strain.png", dpi=400)
plt.close()
print("[+] saved bandgap_SOC_noSOC_vs_strain.png")

print("\n===== ALL DONE =====")
