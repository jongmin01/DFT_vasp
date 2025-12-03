import procar_parser as pp
import numpy as np

# Parse PROCAR
parser = pp.ProcarParser('PROCAR.gz')
data = parser.parse()

print(f"Total k-points: {len(data.kpoints)}")
print(f"\nFirst 10 k-points:")
for i in range(min(10, len(data.kpoints))):
    kpt = data.kpoints[i]
    print(f"  k-point {i:3d}: ({kpt[0]:8.5f}, {kpt[1]:8.5f}, {kpt[2]:8.5f})")

# Look for K and K' points (for hexagonal BZ)
# K point is typically at (1/3, 1/3, 0) in fractional coordinates
# K' point is typically at (2/3, 2/3, 0) or (-1/3, -1/3, 0)

print("\n\nSearching for K and K' points:")
print("K ~ (0.333, 0.333, 0) or similar")
print("K' ~ (0.667, 0.667, 0) or (-0.333, -0.333, 0) or similar\n")

k_candidates = []
kp_candidates = []

for i, kpt in enumerate(data.kpoints):
    # Check for K point (around 1/3, 1/3, 0)
    if abs(kpt[2]) < 0.01:  # z ~ 0
        if abs(abs(kpt[0]) - 1/3) < 0.05 and abs(abs(kpt[1]) - 1/3) < 0.05:
            k_candidates.append((i, kpt))
            print(f"K candidate {i:3d}: ({kpt[0]:8.5f}, {kpt[1]:8.5f}, {kpt[2]:8.5f})")
        
        # Check for K' point (around 2/3, 2/3, 0 or -1/3, -1/3, 0)
        if (abs(abs(kpt[0]) - 2/3) < 0.05 and abs(abs(kpt[1]) - 2/3) < 0.05):
            kp_candidates.append((i, kpt))
            print(f"K' candidate {i:3d}: ({kpt[0]:8.5f}, {kpt[1]:8.5f}, {kpt[2]:8.5f})")

print(f"\n\nLast 10 k-points:")
for i in range(max(0, len(data.kpoints)-10), len(data.kpoints)):
    kpt = data.kpoints[i]
    print(f"  k-point {i:3d}: ({kpt[0]:8.5f}, {kpt[1]:8.5f}, {kpt[2]:8.5f})")
