#!/usr/bin/env bash
# Create SOC-OFF work dirs & inputs for 3-step pipeline per strain dir.
# Steps: 11_scf_nosoc_warm -> 12_scf_nosoc_hr -> 24_band_nosoc_HR
# Run this from the MoS2 root that contains strain_* directories.

set -euo pipefail

echo "[INFO] CWD = $(pwd)"
echo "[INFO] Looking for strain_* directories..."
shopt -s nullglob
STRAINS=(strain_*)
if [ ${#STRAINS[@]} -eq 0 ]; then
  echo "[ERROR] No strain_* directories found in: $(pwd)"
  echo "        Run this script in the folder that has strain_099, strain_100, strain_101, ..."
  exit 1
fi
printf "[INFO] Found strains: %s\n" "${STRAINS[@]}"

# --- config ---
SCF_MESH="24 24 1"     # SCF k-mesh
ENCUT=600
EDIFF_SCF=1E-7
SIGMA=0.05
NEDOS=5000
BAND_N=120             # points per segment for Γ–K–M–Γ

make_kmesh () {
cat <<EOF
Automatic mesh
0
Gamma
$SCF_MESH
0 0 0
EOF
}

make_runpbs () {
cat <<'EOF'
#!/bin/bash -l
#PBS -N VASP_job
#PBS -A cnm84150
#PBS -q batch
#PBS -l nodes=1:ppn=8
#PBS -l walltime=00:40:00
#PBS -j oe
#PBS -V
#PBS -m n

cd "$PBS_O_WORKDIR"
module purge
module load vasp5/5.4/impi-5/intel-16/5.4.1.3-11

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

mpirun -np 8 vasp_std > vasp.out
EOF
}

make_incar_warm () {
cat <<EOF
SYSTEM = MoS2 monolayer SCF-warm (SOC OFF)

# --- Electronic ---
ISTART = 0
ICHARG = 2
ISPIN  = 2
LSORBIT = .FALSE.
LNONCOLLINEAR = .FALSE.

ENCUT  = $ENCUT
EDIFF  = $EDIFF_SCF
ISMEAR = 0
SIGMA  = $SIGMA
PREC   = Accurate
ALGO   = Normal
LREAL  = .FALSE.

# --- Ions/Cell ---
IBRION = -1
NSW    = 0
ISIF   = 2

# --- Slab dipole (2D) ---
LDIPOL = .TRUE.
IDIPOL = 3
DIPOL  = 0.5 0.5 0.5

# --- Outputs ---
LWAVE  = .TRUE.
LCHARG = .TRUE.
NEDOS  = $NEDOS
EOF
}

make_incar_hr () {
cat <<EOF
SYSTEM = MoS2 monolayer SCF-highres (SOC OFF)

# --- Electronic ---
ISTART = 1            # reuse WAVECAR
ICHARG = 11           # fixed charge density (from CHGCAR)
ISPIN  = 2
LSORBIT = .FALSE.
LNONCOLLINEAR = .FALSE.

ENCUT  = $ENCUT
EDIFF  = $EDIFF_SCF
ISMEAR = 0
SIGMA  = $SIGMA
PREC   = Accurate
ALGO   = Normal
LREAL  = .FALSE.
NELM   = 200

# --- Ions/Cell ---
IBRION = -1
NSW    = 0
ISIF   = 2

# --- Slab dipole (2D) ---
LDIPOL = .TRUE.
IDIPOL = 3
DIPOL  = 0.5 0.5 0.5

# --- Outputs ---
LWAVE  = .TRUE.
LCHARG = .TRUE.
NEDOS  = $NEDOS
EOF
}

make_incar_band () {
cat <<EOF
SYSTEM = MoS2 band (SOC OFF)

# --- Electronic (static band) ---
ISTART = 1            # read WAVECAR
ICHARG = 11           # read CHGCAR (fixed)
ISPIN  = 2
LSORBIT = .FALSE.
LNONCOLLINEAR = .FALSE.

ENCUT  = $ENCUT
EDIFF  = 1E-8
ISMEAR = 0
SIGMA  = $SIGMA
PREC   = Accurate
ALGO   = Normal
LREAL  = .FALSE.

# --- No ions/cell movement ---
IBRION = -1
NSW    = 0
ISIF   = 2

# --- Outputs ---
LWAVE  = .FALSE.
LCHARG = .FALSE.
EOF
}

write_band_kpoints () {
cat <<EOF
MoS2 band path (no SOC): Γ–K–M–Γ
$BAND_N
Line-mode
reciprocal
  0.000000 0.000000 0.0   ! Γ
  0.333333 0.333333 0.0   ! K
  0.500000 0.500000 0.0   ! M
  0.000000 0.000000 0.0   ! Γ
EOF
}

# Helper: copy POSCAR/POTCAR from nearest available place
ensure_input_files () {
  local target="$1"
  local src_pos=""
  local src_pot=""

  for cand in "$target/POSCAR" "$(dirname "$target")/POSCAR" "$(dirname "$target")/base/POSCAR" ; do
    [ -f "$cand" ] && { src_pos="$cand"; break; }
  done
  for cand in "$target/POTCAR" "$(dirname "$target")/POTCAR" "$(dirname "$target")/base/POTCAR" ; do
    [ -f "$cand" ] && { src_pot="$cand"; break; }
  done

  if [ -n "$src_pos" ]; then
    cp -f "$src_pos" "$target/POSCAR"
  else
    echo "[WARN] POSCAR not found for $target (looked in self/..//../base)."
  fi
  if [ -n "$src_pot" ]; then
    cp -f "$src_pot" "$target/POTCAR"
  else
    echo "[WARN] POTCAR not found for $target (looked in self/..//../base)."
  fi
}

# --- per-strain setup ---
for d in "${STRAINS[@]}"; do
  [ -d "$d" ] || continue
  echo "[MAKE] $d"

  # step1
  mkdir -p "$d/11_scf_nosoc_warm"
  ( cd "$d/11_scf_nosoc_warm"
    make_incar_warm > INCAR
    make_kmesh     > KPOINTS
    make_runpbs    > run.pbs
    chmod +x run.pbs
    ensure_input_files "$(pwd)"
  )

  # step2
  mkdir -p "$d/12_scf_nosoc_hr"
  ( cd "$d/12_scf_nosoc_hr"
    make_incar_hr > INCAR
    make_kmesh    > KPOINTS
    make_runpbs   > run.pbs
    chmod +x run.pbs
    ensure_input_files "$(pwd)"
  )

  # step3
  mkdir -p "$d/24_band_nosoc_HR"
  ( cd "$d/24_band_nosoc_HR"
    make_incar_band   > INCAR
    write_band_kpoints > KPOINTS
    make_runpbs       > run.pbs
    chmod +x run.pbs
    ensure_input_files "$(pwd)"
  )
done

echo "[DONE] SOC-OFF workdirs ready."
