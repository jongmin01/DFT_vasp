import numpy as np
import procar_parser as pp
import plot_band_pdos as pbp
import plot_layer_resolved as plr
from scipy.ndimage import gaussian_filter1d

# 1. 파일 설정
procar_file = "PROCAR.gz"
poscar_file = "CONTCAR"  # 또는 POSCAR
doscar_file = "DOSCAR"
output_file = "final_layer_band_pdos.png"

print("1. 데이터 파싱 중...")
parser = pp.ProcarParser(procar_file)
data = parser.parse()

print("2. 층 정보 인식 중...")
layer_atoms = plr.identify_tmdc_layers(poscar_file)
target_layers = ['MoSSe', 'WSSe']
# 실제 원자가 있는 층만 선택
layer_groups = {k: layer_atoms[k] for k in target_layers if k in layer_atoms}

# 3. 밴드 프로젝션 (왼쪽 패널용)
print("3. 밴드 프로젝션 계산 중...")
layer_projections = {}
for name, atoms in layer_groups.items():
    proj = parser.get_atom_projection(atoms)
    if data.is_spin_polarized:
        proj = proj[:, :, 0]  # Spin-up만 사용 (nkpts, nbands)
    layer_projections[name] = proj

# 4. PDOS 계산 (오류 수정됨: 직접 계산 방식)
print("4. PDOS 계산 중 (Direct Method)...")

# 에너지 그리드 생성 (밴드 에너지 범위 기반)
energies = data.energies[:, :, 0] if data.is_spin_polarized else data.energies
e_min, e_max = energies.min(), energies.max()
e_grid = np.linspace(e_min, e_max, 2000) # 촘촘하게
e_centers = (e_grid[:-1] + e_grid[1:]) / 2

# DOSCAR에서 페르미 레벨 가져오기
_, _, efermi = pbp.read_doscar_total(doscar_file)

pdos_dict = {}
for name, atoms in layer_groups.items():
    # 1) 안전한 parser 메서드로 projection 가져오기
    proj = parser.get_atom_projection(atoms)
    
    # 2) Spin-up 성분만 추출
    if data.is_spin_polarized:
        proj = proj[:, :, 0]
    
    # 3) 히스토그램으로 DOS 계산 (이제 shape이 (nkpts, nbands)로 확실히 일치함)
    dos, _ = np.histogram(energies.flatten(), bins=e_grid, 
                          weights=proj.flatten(), density=True)
    
    # 4) 부드럽게 만들기 (Gaussian Smearing)
    dos_smooth = gaussian_filter1d(dos, sigma=5) # sigma 조절로 부드러움 변경 가능
    pdos_dict[name] = dos_smooth

# 5. 그림 그리기
print(f"5. 그림 저장 중... ({output_file})")

# plot_band_pdos.py의 그리기 함수 호출
pbp.plot_band_with_layer_pdos(
    kpoints=data.kpoints,
    energies=energies,
    layer_projections=layer_projections,
    dos_energy=e_centers,
    pdos_dict=pdos_dict,
    efermi=efermi,
    title="MoSSe/WSSe Heterostructure Band & PDOS",
    output=output_file,
    energy_range=(-3, 3)
)

print("완료! 확인해보세요.")
