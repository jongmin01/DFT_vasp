#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
build_correct_heterostructure.py
---------------------------------
올바른 Janus TMDC heterostructure 생성
POTCAR 순서: Mo W S Se
"""

import os
import numpy as np
from ase.io import read, write
from ase import Atoms

# ============================================================
# 설정
# ============================================================

# 입력 파일 (seed 구조)
input_files = {
    'MoSSe_S_up': '/mnt/user-data/uploads/MoSeS_S_up.vasp',
    'MoSSe_Se_up': '/mnt/user-data/uploads/MoSeS_Se_up.vasp',
    'WSSe_S_up': '/mnt/user-data/uploads/WSeS_S_up.vasp',
    'WSSe_Se_up': '/mnt/user-data/uploads/WSeS_Se_up.vasp'
}

# Heterostructure 조합
combinations = [
    ('MoSSe_S_up', 'WSSe_S_up'),
    ('MoSSe_S_up', 'WSSe_Se_up'),
    ('MoSSe_Se_up', 'WSSe_S_up'),
    ('MoSSe_Se_up', 'WSSe_Se_up'),
]

# 구조 파라미터
supercell = (5, 5, 1)           # Supercell 크기
interlayer_gap = 6.8            # 층간 거리 (Å)
vacuum = 15.0                   # Vacuum 공간 (Å)
potcar_order = ['Mo', 'W', 'S', 'Se']  # POTCAR 순서

output_dir = '/mnt/user-data/outputs/correct_heterostructures'

# ============================================================
# 함수들
# ============================================================

def fix_janus_structure(atoms):
    """
    Periodic boundary를 고려하여 Janus 구조를 올바르게 재배치
    X2-M-X1 구조를 컴팩트하게 만듦 (X2가 periodic으로 아래에 있을 경우 처리)
    """
    positions = atoms.get_positions()
    symbols = atoms.get_chemical_symbols()
    cell = atoms.get_cell()
    c_axis = cell[2, 2]
    
    # 원자 분류
    metal_idx = None
    chalcogen_data = []  # (index, symbol, z_position)
    
    for i, sym in enumerate(symbols):
        if sym in ['Mo', 'W']:
            metal_idx = i
        elif sym in ['S', 'Se']:
            chalcogen_data.append((i, sym, positions[i, 2]))
    
    if metal_idx is None or len(chalcogen_data) != 2:
        raise ValueError(f"Invalid Janus structure: need 1 metal + 2 chalcogens, got {symbols}")
    
    metal_sym = symbols[metal_idx]
    metal_z = positions[metal_idx, 2]
    
    # Chalcogen들을 periodic wrapping 고려하여 재배치
    wrapped_chalcogens = []
    for idx, sym, z in chalcogen_data:
        # z > c/2인 원자는 periodic으로 아래쪽에 있을 가능성
        if z > c_axis / 2:
            z_wrapped = z - c_axis
        else:
            z_wrapped = z
        wrapped_chalcogens.append((idx, sym, z_wrapped))
    
    # Metal 기준 상대 위치로 재배치
    new_positions = np.zeros((3, 3))
    new_symbols = []
    
    # Metal을 중간(z=1.6)에 배치
    metal_relative_z = 1.6
    new_positions[metal_idx] = [
        positions[metal_idx, 0],
        positions[metal_idx, 1],
        metal_relative_z
    ]
    new_symbols.append(metal_sym)
    
    # Chalcogen들을 metal 기준으로 배치
    chal_positions = []
    for idx, sym, z_wrapped in wrapped_chalcogens:
        rel_dist = z_wrapped - metal_z  # Metal 기준 상대 거리
        new_z = metal_relative_z + rel_dist
        
        new_positions[idx] = [
            positions[idx, 0],
            positions[idx, 1],
            new_z
        ]
        chal_positions.append((idx, sym, new_z, rel_dist))
    
    # 모든 원자를 양수 z로 shift
    z_coords = new_positions[:, 2]
    z_min = z_coords.min()
    new_positions[:, 2] -= z_min
    
    # 층 두께
    thickness = new_positions[:, 2].max()
    
    # 새로운 cell (compact + 10Å vacuum)
    new_cell = cell.copy()
    new_cell[2, 2] = thickness + 10.0
    
    new_atoms = Atoms(
        symbols=symbols,
        positions=new_positions,
        cell=new_cell,
        pbc=True
    )
    
    return new_atoms


def build_heterostructure(bot_name, top_name):
    """두 층으로 heterostructure 생성"""
    
    print(f"\n{'='*70}")
    print(f"Building: {bot_name} // {top_name}")
    print(f"{'='*70}")
    
    # 입력 파일 읽기
    bot_atoms = read(input_files[bot_name])
    top_atoms = read(input_files[top_name])
    
    print(f"  원본 구조:")
    print(f"    Bottom: {len(bot_atoms)} atoms, c={bot_atoms.cell[2,2]:.4f} Å")
    print(f"    Top: {len(top_atoms)} atoms, c={top_atoms.cell[2,2]:.4f} Å")
    
    # Janus 구조 수정
    print(f"  Janus 구조 수정 중...")
    bot_fixed = fix_janus_structure(bot_atoms)
    top_fixed = fix_janus_structure(top_atoms)
    
    bot_thickness = bot_fixed.positions[:, 2].max() - bot_fixed.positions[:, 2].min()
    top_thickness = top_fixed.positions[:, 2].max() - top_fixed.positions[:, 2].min()
    
    print(f"    Bottom fixed: thickness={bot_thickness:.4f} Å")
    print(f"    Top fixed: thickness={top_thickness:.4f} Å")
    
    # Supercell 생성
    nx, ny, nz = supercell
    bot_super = bot_fixed.repeat((nx, ny, nz))
    top_super = top_fixed.repeat((nx, ny, nz))
    
    print(f"  Supercell {nx}x{ny}x{nz}:")
    print(f"    Bottom: {len(bot_super)} atoms")
    print(f"    Top: {len(top_super)} atoms")
    
    # 위치 조정
    bot_positions = bot_super.get_positions()
    top_positions = top_super.get_positions()
    
    # Bottom layer를 z=0부터 시작
    bot_positions[:, 2] -= bot_positions[:, 2].min()
    bot_z_max = bot_positions[:, 2].max()
    
    # Top layer를 gap만큼 띄워서 배치
    top_positions[:, 2] -= top_positions[:, 2].min()
    top_positions[:, 2] += (bot_z_max + interlayer_gap)
    top_z_max = top_positions[:, 2].max()
    
    print(f"  층 배치:")
    print(f"    Bottom: 0.000 → {bot_z_max:.4f} Å")
    print(f"    Gap: {interlayer_gap:.4f} Å")
    print(f"    Top: {bot_z_max + interlayer_gap:.4f} → {top_z_max:.4f} Å")
    
    # 전체 cell 높이
    total_height = top_z_max + vacuum
    
    # 원자 결합
    all_symbols = bot_super.get_chemical_symbols() + top_super.get_chemical_symbols()
    all_positions = np.vstack([bot_positions, top_positions])
    
    # 새 cell 생성
    new_cell = bot_super.cell.copy()
    new_cell[2, 2] = total_height
    
    hetero = Atoms(
        symbols=all_symbols,
        positions=all_positions,
        cell=new_cell,
        pbc=True
    )
    
    print(f"  전체 구조:")
    print(f"    총 원자 수: {len(hetero)}")
    print(f"    Cell c: {total_height:.4f} Å")
    print(f"    Vacuum: {vacuum:.4f} Å")
    
    return hetero


def write_poscar_with_potcar_order(filename, atoms, comment, potcar_order):
    """
    POTCAR 순서에 맞게 POSCAR 작성
    potcar_order = ['Mo', 'W', 'S', 'Se']
    """
    symbols = np.array(atoms.get_chemical_symbols())
    positions = atoms.get_scaled_positions()
    cell = atoms.get_cell()
    
    # 원소별로 분류
    sorted_indices = []
    element_counts = []
    present_elements = []
    
    for elem in potcar_order:
        indices = np.where(symbols == elem)[0]
        if len(indices) > 0:
            sorted_indices.extend(indices)
            element_counts.append(len(indices))
            present_elements.append(elem)
    
    # 정렬된 순서대로 재배치
    sorted_positions = positions[sorted_indices]
    sorted_symbols = symbols[sorted_indices]
    
    # POSCAR 작성
    with open(filename, 'w') as f:
        f.write(f"{comment}\n")
        f.write("1.0\n")
        
        # Lattice vectors
        for vec in cell:
            f.write(f"  {vec[0]:20.16f}  {vec[1]:20.16f}  {vec[2]:20.16f}\n")
        
        # Element names (POTCAR 순서)
        f.write("  " + "  ".join(present_elements) + "\n")
        
        # Element counts
        f.write("  " + "  ".join(map(str, element_counts)) + "\n")
        
        # Direct coordinates
        f.write("Direct\n")
        for pos in sorted_positions:
            f.write(f"  {pos[0]:20.16f}  {pos[1]:20.16f}  {pos[2]:20.16f}\n")
    
    return present_elements, element_counts


def verify_structure(atoms, name):
    """구조 검증"""
    print(f"\n  {'='*60}")
    print(f"  구조 검증: {name}")
    print(f"  {'='*60}")
    
    positions = atoms.get_positions()
    symbols = atoms.get_chemical_symbols()
    
    # z 좌표로 정렬하여 층 구조 확인
    z_order = np.argsort(positions[:, 2])
    
    # 고유한 z 레벨 찾기 (tolerance: 0.1 Å)
    z_levels = []
    current_level = []
    current_z = positions[z_order[0], 2]
    
    for idx in z_order:
        z = positions[idx, 2]
        if abs(z - current_z) < 0.1:
            current_level.append(symbols[idx])
        else:
            z_levels.append((current_z, current_level))
            current_level = [symbols[idx]]
            current_z = z
    z_levels.append((current_z, current_level))
    
    print(f"  층 구조 (아래 → 위):")
    for i, (z, level_atoms) in enumerate(z_levels, 1):
        from collections import Counter
        counts = Counter(level_atoms)
        atom_str = ", ".join([f"{elem}×{cnt}" for elem, cnt in sorted(counts.items())])
        print(f"    {i}. z={z:7.4f} Å: {atom_str}")
    
    # Metal 원자들의 결합 확인
    for metal in ['Mo', 'W']:
        metal_indices = [i for i, sym in enumerate(symbols) if sym == metal]
        if not metal_indices:
            continue
        
        # 첫 번째 metal 원자로 대표 확인
        metal_idx = metal_indices[0]
        metal_pos = positions[metal_idx]
        
        # 가장 가까운 chalcogen 2개 찾기
        distances = []
        for i, sym in enumerate(symbols):
            if sym in ['S', 'Se']:
                dist = np.linalg.norm(metal_pos - positions[i])
                distances.append((dist, sym, positions[i, 2]))
        
        distances.sort()
        if len(distances) >= 2:
            print(f"\n  {metal} 원자 결합:")
            print(f"    1. {distances[0][1]:3s}: {distances[0][0]:.4f} Å (z={distances[0][2]:.4f})")
            print(f"    2. {distances[1][1]:3s}: {distances[1][0]:.4f} Å (z={distances[1][2]:.4f})")
            
            if distances[0][0] < 3.0 and distances[1][0] < 3.0:
                if distances[0][1] != distances[1][1]:
                    print(f"    ✓ Janus 구조 확인!")
                else:
                    print(f"    ⚠ 동일한 chalcogen과만 결합")
            else:
                print(f"    ✗ 결합 거리 이상")


# ============================================================
# Main
# ============================================================

def main():
    """Main execution"""
    
    print("\n" + "="*70)
    print("올바른 Janus TMDC Heterostructure 생성기")
    print("="*70)
    print(f"POTCAR 순서: {' '.join(potcar_order)}")
    print(f"Supercell: {supercell[0]}×{supercell[1]}×{supercell[2]}")
    print(f"층간 거리: {interlayer_gap} Å")
    print(f"Vacuum: {vacuum} Å")
    print("="*70)
    
    # 출력 디렉토리 생성
    os.makedirs(output_dir, exist_ok=True)
    
    # 각 조합별로 생성
    success_count = 0
    for bot_name, top_name in combinations:
        try:
            # Heterostructure 생성
            hetero = build_heterostructure(bot_name, top_name)
            
            # 구조 검증
            verify_structure(hetero, f"{bot_name}//{top_name}")
            
            # 출력 디렉토리
            dir_name = f"{bot_name}__{top_name}"
            output_path = os.path.join(output_dir, dir_name)
            os.makedirs(output_path, exist_ok=True)
            
            # POSCAR 작성 (POTCAR 순서)
            poscar_file = os.path.join(output_path, 'POSCAR')
            comment = f"{bot_name}__{top_name} heterostructure"
            elements, counts = write_poscar_with_potcar_order(
                poscar_file, hetero, comment, potcar_order
            )
            
            print(f"\n  ✓ POSCAR 생성: {poscar_file}")
            print(f"    원소 순서: {' '.join(elements)}")
            print(f"    원소 개수: {' '.join(map(str, counts))}")
            print(f"    POTCAR 순서와 일치: ✓")
            
            success_count += 1
            
        except Exception as e:
            print(f"\n  ✗ 오류 발생: {e}")
            import traceback
            traceback.print_exc()
    
    # 최종 요약
    print("\n" + "="*70)
    print(f"완료: {success_count}/{len(combinations)} 구조 생성")
    print("="*70)
    print(f"\n출력 디렉토리: {output_dir}")
    print("\n다음 단계:")
    print("  1. POSCAR 검증:")
    print(f"     head -10 {output_dir}/*/POSCAR")
    print("  2. 시각화 확인")
    print("  3. Carbon으로 전송")
    print("="*70 + "\n")


if __name__ == '__main__':
    main()
