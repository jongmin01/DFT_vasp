#!/usr/bin/env bash
# Create ±1% in-plane strain sets and populate SOC SCF/Band inputs.
# Requires: base/POSCAR, base/POTCAR present.

set -euo pipefail

BASE=base
[[ -f $BASE/POSCAR && -f $BASE/POTCAR ]] || { echo "Missing base/POSCAR or base/POTCAR"; exit 1; }

# Strain factors (scale a,b only; keep c)
STRAINS=("0.99" "1.00" "1.01")
MAPIDX=("099" "100" "101")  # for nicer directory names

for i in "${!STRAINS[@]}"; do
  s=${STRAINS[$i]}
  tag=${MAPIDX[$i]}
  root="strain_${tag}"
  echo ">> Making ${root} (scale a,b by ${s})"
  rm -rf "${root}"; mkdir -p "${root}"

  # --- Build POSCAR for this strain: scale first two lattice vectors only ---
  awk -v s="$s" 'NR==1{print;next}
    NR==2{print;next}
    NR==3{printf("%20.16f %20.16f %20.16f\n",$1*s,$2*s,$3*s); next}
    NR==4{printf("%20.16f %20.16f %20.16f\n",$1*s,$2*s,$3*s); next}
    {print}' $BASE/POSCAR > "${root}/POSCAR"

  cp $BASE/POTCAR "${root}/POTCAR"

  # --- Stage directories ---
  mkdir -p "${root}/11_scf_soc_warm" "${root}/12_scf_soc_hr" "${root}/24_band_soc_HR"

  # Common files
  for d in 11_scf_soc_warm 12_scf_soc_hr 24_band_soc_HR; do
    cp "${root}/POSCAR" "${root}/${d}/POSCAR"
    cp "${root}/POTCAR" "${root}/${d}/POTCAR"
  done

  # -------- INCAR templates --------
  # SCF step-1 (warm, ISTART=0, 12x12x1)
  cat > "${root}/11_scf_soc_warm/INCAR" <<'EOF'
SYSTEM = MoS2 SOC SCF (warm, 12x12x1)
ISTART = 0
ICHARG = 11
ENCUT  = 600
EDIFF  = 1E-7
ISMEAR = 0
SIGMA  = 0.05
PREC   = Accurate
ALGO   = All
LREAL  = .FALSE.
LASPH  = .TRUE.
NELM   = 200
LMAXMIX  = 4
AMIX     = 0.2
BMIX     = 0.0001
AMIX_MAG = 0.2
BMIX_MAG = 0.0001
LSORBIT  = .TRUE.
GGA_COMPAT = .FALSE.
ISYM     = 0
SAXIS    = 0 0 1
IBRION = -1
NSW    = 0
ISIF   = 2
LDIPOL = .TRUE.
IDIPOL = 3
DIPOL  = 0.5 0.5 0.5
LWAVE  = .TRUE.
LCHARG = .TRUE.
NWRITE = 3
EOF

  cat > "${root}/11_scf_soc_warm/KPOINTS" <<'EOF'
Automatic mesh
0
Gamma
12 12 1
0 0 0
EOF

  # SCF step-2 (HR, ISTART=1, 24x24x1)
  cat > "${root}/12_scf_soc_hr/INCAR" <<'EOF'
SYSTEM = MoS2 SOC SCF (high-res, 24x24x1)
ISTART = 1
ICHARG = 11
ENCUT  = 600
EDIFF  = 1E-7
ISMEAR = 0
SIGMA  = 0.03
PREC   = Accurate
ALGO   = All
LREAL  = .FALSE.
LASPH  = .TRUE.
NELM   = 200
LMAXMIX  = 4
AMIX     = 0.2
BMIX     = 0.0001
AMIX_MAG = 0.2
BMIX_MAG = 0.0001
LSORBIT  = .TRUE.
GGA_COMPAT = .FALSE.
ISYM     = 0
SAXIS    = 0 0 1
IBRION = -1
NSW    = 0
ISIF   = 2
LDIPOL = .TRUE.
IDIPOL = 3
DIPOL  = 0.5 0.5 0.5
LWAVE  = .TRUE.
LCHARG = .TRUE.
NWRITE = 3
EOF

  cat > "${root}/12_scf_soc_hr/KPOINTS" <<'EOF'
Automatic mesh
0
Gamma
24 24 1
0 0 0
EOF

  # Band (SOC, line-mode N=120, NBANDS ample)
  cat > "${root}/24_band_soc_HR/INCAR" <<'EOF'
SYSTEM = MoS2 Band + SOC (HR, N=120)
ISTART = 1
ICHARG = 11
ENCUT  = 600
EDIFF  = 1E-8
ISMEAR = 0
SIGMA  = 0.05
PREC   = Accurate
ALGO   = Normal
LREAL  = .FALSE.
LASPH  = .TRUE.
LSORBIT  = .TRUE.
GGA_COMPAT = .FALSE.
ISYM     = 0
SAXIS    = 0 0 1
IBRION = -1
NSW    = 0
ISIF   = 2
LDIPOL = .TRUE.
IDIPOL = 3
DIPOL  = 0.5 0.5 0.5
LWAVE  = .FALSE.
LCHARG = .FALSE.
NWRITE = 3
NBANDS = 96
EOF

  cat > "${root}/24_band_soc_HR/KPOINTS" <<'EOF'
MoS2 band path (SOC): Γ–K–M–Γ (120 pts/seg)
120
Line-mode
reciprocal
 0.000000 0.000000 0.0   ! Γ
 0.333333 0.333333 0.0   ! K

 0.333333 0.333333 0.0   ! K
 0.500000 0.500000 0.0   ! M

 0.500000 0.500000 0.0   ! M
 0.000000 0.000000 0.0   ! Γ
EOF

done

echo "All strain workdirs prepared: ${STRAINS[*]}"
