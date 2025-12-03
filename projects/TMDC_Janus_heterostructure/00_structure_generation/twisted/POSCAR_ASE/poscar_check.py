import numpy as np

with open('POSCAR', 'r') as f:
    lines = f.readlines()

# Read lattice
scale = float(lines[1].strip())
lattice = []
for i in range(2, 5):
    vec = [float(x) for x in lines[i].split()]
    lattice.append(vec)
lattice = np.array(lattice) * scale

# Read positions
elements = lines[5].split()
counts = [int(x) for x in lines[6].split()]
total = sum(counts)

positions = []
for i in range(8, 8 + total):
    pos = [float(x) for x in lines[i].split()[:3]]
    positions.append(pos)

positions = np.array(positions)
positions_cart = positions @ lattice

z_coords = positions_cart[:, 2]
z_sorted = np.sort(z_coords)

# Find layer gap
z_gaps = np.diff(z_sorted)
max_gap_idx = np.argmax(z_gaps)

layer1_max = z_sorted[max_gap_idx]
layer2_min = z_sorted[max_gap_idx + 1]
interlayer = layer2_min - layer1_max

print(f"Interlayer distance: {interlayer:.3f} Å")
print(f"Layer 1 z-range: {z_sorted[0]:.3f} - {layer1_max:.3f} Å")
print(f"Layer 2 z-range: {layer2_min:.3f} - {z_sorted[-1]:.3f} Å")

# Check for atomic overlaps
min_distance = 100
for i in range(len(positions_cart)):
    for j in range(i+1, len(positions_cart)):
        dist = np.linalg.norm(positions_cart[i] - positions_cart[j])
        if dist < min_distance:
            min_distance = dist

print(f"\nMinimum atomic distance: {min_distance:.3f} Å")
if min_distance < 1.5:
    print("⚠️  WARNING: Atoms too close! Possible overlap!")
