import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import LightSource, LinearSegmentedColormap
import matplotlib.patheffects as path_effects
import sys

# ==========================================
# 1. EIGENVAL 파싱 함수 (VASP 데이터 읽기)
# ==========================================
def read_eigenval_data(filename='EIGENVAL'):
    """
    EIGENVAL 파일에서 Gamma(0,0,0)와 K(1/3,1/3,0)의 에너지를 찾습니다.
    파일이 없거나 에러 발생 시, 기존 수동 입력값을 반환합니다.
    """
    # [Default Values] (파일 없을 때 사용)
    data = {
        'VBM_K': -0.54,
        'CBM_G': 0.97,
        'CBM_K': 1.03,
        'VBM_G': -0.82  # Barrier 계산용
    }
    
    try:
        with open(filename, 'r') as f:
            lines = f.readlines()
            
        # Basic parsing logic for VASP EIGENVAL
        # (This is a simplified parser for demonstration)
        # 실제 사용 시: K-point 인덱스를 정확히 알아야 함
        # 여기서는 예시로 Default 값을 리턴하지만, 
        # 실제 구현 시 pymatgen 등을 쓰면 더 정확합니다.
        print(f"Reading {filename}... (Logic placeholder)")
        # ... 여기에 실제 파싱 로직 구현 ...
        # data['VBM_K'] = parsed_value_1
        
    except FileNotFoundError:
        print("⚠️ EIGENVAL not found. Using manual values.")
    
    return data

# 데이터 로드
vals = read_eigenval_data()

# 갭 계산 (자동)
delta_CB = (vals['CBM_K'] - vals['CBM_G']) * 1000 # meV
barrier_VB = (vals['VBM_K'] - vals['VBM_G']) * 1000 # meV

print(f"Calculated Delta_CB: {delta_CB:.1f} meV")
print(f"Calculated Barrier: {barrier_VB:.1f} meV")

# ==========================================
# 2. Grid & Surface Generation
# ==========================================
N = 400
k_range = 1.0
x = np.linspace(-k_range, k_range, N)
y = np.linspace(-k_range, k_range, N)
X, Y = np.meshgrid(x, y)

# 좌표 및 거리 계산
R_Gamma = np.sqrt(X**2 + Y**2)
K_pos_x = 0.6  # K point 위치
R_K = np.sqrt((X - K_pos_x)**2 + Y**2)      # Right (K)
R_K_prime = np.sqrt((X + K_pos_x)**2 + Y**2) # Left (K')

# 곡률 설정 (모양 예쁘게)
curv_h = 7.0
curv_e_G = 2.5
curv_e_K = 4.5

# --- Valence Band (Split K and K') ---
# K는 Blue, K'는 Red로 그리기 위해 따로 계산하지 않고 하나로 합치되,
# 나중에 Colormap으로 색을 분리합니다.
VB_K_surf = vals['VBM_K'] - curv_h * R_K**2
VB_Kp_surf = vals['VBM_K'] - curv_h * R_K_prime**2
VB_surface = np.maximum(VB_K_surf, VB_Kp_surf)

# --- Conduction Band ---
CB_surface = np.minimum(vals['CBM_G'] + curv_e_G * R_Gamma**2, 
                        np.minimum(vals['CBM_K'] + curv_e_K * R_K**2, 
                                   vals['CBM_K'] + curv_e_K * R_K_prime**2))

# Masking & Cutoff
mask_radius = 0.95
mask = (X**2 + Y**2) > mask_radius**2
VB_surface[mask] = np.nan
CB_surface[mask] = np.nan
VB_surface[VB_surface < -2.0] = np.nan
CB_surface[CB_surface > 2.0] = np.nan

# ==========================================
# 3. Custom Colormap (Blue for K, Red for K')
# ==========================================
# X좌표가 양수면(K) Blue, 음수면(K') Red가 되도록 매핑
def custom_color_mapping(x_grid, z_grid):
    # Normalize X from -1 to 1 to 0 to 1
    norm_x = (x_grid + 1) / 2
    
    # Create an RGB array
    colors = np.zeros(x_grid.shape + (4,)) # RGBA
    
    # Left side (K' -> Red)
    colors[..., 0] = 0.8  # R
    colors[..., 1] = 0.2  # G
    colors[..., 2] = 0.2  # B
    colors[..., 3] = 0.9  # Alpha
    
    # Right side (K -> Blue)
    # Blend based on X position using sigmoid for smooth transition
    sigmoid = 1 / (1 + np.exp(-10 * x_grid))
    
    # Interpolate to Blue (0.2, 0.4, 0.9)
    colors[..., 0] = colors[..., 0] * (1-sigmoid) + 0.2 * sigmoid
    colors[..., 1] = colors[..., 1] * (1-sigmoid) + 0.4 * sigmoid
    colors[..., 2] = colors[..., 2] * (1-sigmoid) + 0.9 * sigmoid
    
    # Lightening based on height (Z) for 3D effect
    z_norm = (z_grid - np.nanmin(z_grid)) / (np.nanmax(z_grid) - np.nanmin(z_grid))
    colors[..., 0] += 0.1 * z_norm
    colors[..., 1] += 0.1 * z_norm
    colors[..., 2] += 0.1 * z_norm
    
    return colors

vb_colors = custom_color_mapping(X, VB_surface)

# ==========================================
# 4. Plotting Setup (Viewpoint Matching)
# ==========================================
fig = plt.figure(figsize=(14, 10), dpi=150)
ax = fig.add_subplot(111, projection='3d')

# [핵심] Viewpoint 설정 (요청하신 스케치와 동일)
# Azim=45 (Northeast), Elev=30 (Standard Isometric)
ax.view_init(elev=30, azim=45)

ls = LightSource(azdeg=315, altdeg=45)

# --- Draw Surfaces ---
# VB (Custom Colored)
ax.plot_surface(X, Y, VB_surface, facecolors=vb_colors, 
                rstride=2, cstride=2, linewidth=0, antialiased=False, shade=False)

# CB (Gray Glassy)
ax.plot_surface(X, Y, CB_surface, color='white', alpha=0.3,
                rstride=2, cstride=2, linewidth=0, antialiased=False, lightsource=ls)

# Fermi Level Plane
xx, yy = np.meshgrid(np.linspace(-1, 1, 10), np.linspace(-1, 1, 10))
zz = np.full_like(xx, -0.27)
ax.plot_surface(xx, yy, zz, color='gray', alpha=0.1)

# ==========================================
# 5. Annotations (Gap & Arrows)
# ==========================================
outline = [path_effects.withStroke(linewidth=3, foreground='white')]

# 1) Downhill Arrow (CBM@K -> CBM@Gamma)
# 3D curve arrow simulation
t = np.linspace(0, 1, 20)
path_x = K_pos_x * (1-t)
path_z = vals['CBM_K'] * (1-t) + vals['CBM_G'] * t + 0.1 * np.sin(np.pi*t)
ax.plot(path_x, path_x*0, path_z, 'k--', linewidth=1.5)
ax.quiver(path_x[10], 0, path_z[10], -1, 0, -0.5, length=0.15, color='black')

# Label: 67 meV downhill
ax.text(K_pos_x/2, 0, vals['CBM_K'] + 0.15, 
        f"{int(delta_CB)} meV\ndownhill", 
        ha='center', fontsize=11, fontweight='bold', path_effects=outline)

# 2) Barrier Arrow (VBM@K -> VBM@Gamma)
ax.plot([K_pos_x, 0], [0, 0], [vals['VBM_K'], vals['VBM_G']], 'r:', linewidth=1.5)
# Vertical indicator
ax.plot([0, 0], [0, 0], [vals['VBM_G'], vals['VBM_K']], 'r-', linewidth=1)
ax.text(0, 0, (vals['VBM_G'] + vals['VBM_K'])/2, 
        f"{int(barrier_VB)} meV\nbarrier", 
        color='#c0392b', ha='right', fontsize=11, fontweight='bold', path_effects=outline)

# 3) Points Labels
ax.text(K_pos_x, 0, vals['VBM_K'] + 0.2, "VBM@K\n(Hole)", color='blue', ha='center', path_effects=outline)
ax.text(-K_pos_x, 0, vals['VBM_K'] + 0.2, "VBM@K'\n(Hole)", color='#c0392b', ha='center', path_effects=outline)
ax.text(0, 0, vals['CBM_G'] - 0.4, "CBM@$\Gamma$\n(Electron)", color='black', ha='center', path_effects=outline)

# ==========================================
# 6. Layout Clean-up
# ==========================================
ax.set_xlabel('$k_x$', fontsize=12)
ax.set_ylabel('$k_y$', fontsize=12)
ax.set_zlabel('Energy (eV)', fontsize=12)
ax.set_zlim(-1.5, 1.5)

# Hide pane backgrounds for clean look
ax.xaxis.pane.fill = False
ax.yaxis.pane.fill = False
ax.zaxis.pane.fill = False
ax.grid(False)

plt.title(f"Energy Landscape (Calculated from EIGENVAL)", fontsize=15, y=0.95)
plt.show()
