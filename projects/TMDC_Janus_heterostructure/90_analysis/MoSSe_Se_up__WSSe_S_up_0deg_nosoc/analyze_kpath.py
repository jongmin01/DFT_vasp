import procar_parser as pp
import numpy as np

# Parse PROCAR
parser = pp.ProcarParser('PROCAR.gz')
data = parser.parse()

print("=" * 70)
print("K-PATH ANALYSIS")
print("=" * 70)

# Identify path segments
print("\nAnalyzing k-path structure...")
print(f"Total k-points: {len(data.kpoints)}")

# Find discontinuities (path segments)
segments = []
current_segment = [0]

for i in range(1, len(data.kpoints)):
    dk = np.linalg.norm(data.kpoints[i] - data.kpoints[i-1])
    if dk > 0.15:  # Large jump indicates new segment
        segments.append(current_segment)
        current_segment = [i]
    else:
        current_segment.append(i)

segments.append(current_segment)

print(f"\nFound {len(segments)} path segments:")
for idx, seg in enumerate(segments):
    kstart = data.kpoints[seg[0]]
    kend = data.kpoints[seg[-1]]
    print(f"\nSegment {idx+1}: k-points {seg[0]}-{seg[-1]} ({len(seg)} points)")
    print(f"  Start: ({kstart[0]:7.4f}, {kstart[1]:7.4f}, {kstart[2]:7.4f})")
    print(f"  End:   ({kend[0]:7.4f}, {kend[1]:7.4f}, {kend[2]:7.4f})")
    
    # Check if this segment contains K or K' valley
    for kpt_idx in seg:
        kpt = data.kpoints[kpt_idx]
        # K valley (1/3, 1/3, 0)
        if abs(kpt[2]) < 0.01:
            if abs(kpt[0] - 1/3) < 0.02 and abs(kpt[1] - 1/3) < 0.02:
                print(f"  >>> Contains K valley at k-point {kpt_idx}: ({kpt[0]:.5f}, {kpt[1]:.5f}, {kpt[2]:.5f})")
            # K' valley could be at (2/3, 2/3, 0) or (-1/3, -1/3, 0)
            if (abs(kpt[0] - 2/3) < 0.02 and abs(kpt[1] - 2/3) < 0.02) or \
               (abs(kpt[0] + 1/3) < 0.02 and abs(kpt[1] + 1/3) < 0.02):
                print(f"  >>> Contains K' valley at k-point {kpt_idx}: ({kpt[0]:.5f}, {kpt[1]:.5f}, {kpt[2]:.5f})")

# Common high-symmetry points in hexagonal BZ
print("\n" + "=" * 70)
print("CHECKING FOR STANDARD HIGH-SYMMETRY POINTS")
print("=" * 70)
print("\nStandard points in hexagonal BZ:")
print("  Γ (Gamma): (0, 0, 0)")
print("  M: (0.5, 0, 0) or (0, 0.5, 0)")
print("  K: (1/3, 1/3, 0) or (2/3, 1/3, 0)")
print("  K': (2/3, 2/3, 0) or (-1/3, -1/3, 0) [time-reversal of K]")

found_points = {}
for i, kpt in enumerate(data.kpoints):
    if np.linalg.norm(kpt) < 0.01:
        found_points.setdefault('Γ', []).append(i)
    if abs(kpt[2]) < 0.01:
        if abs(kpt[0] - 0.5) < 0.02 or abs(kpt[1] - 0.5) < 0.02:
            found_points.setdefault('M', []).append(i)
        if abs(kpt[0] - 1/3) < 0.02 and abs(kpt[1] - 1/3) < 0.02:
            found_points.setdefault('K', []).append(i)
        if (abs(kpt[0] - 2/3) < 0.02 and abs(kpt[1] - 2/3) < 0.02) or \
           (abs(kpt[0] + 1/3) < 0.02 and abs(kpt[1] + 1/3) < 0.02):
            found_points.setdefault("K'", []).append(i)

print("\nFound high-symmetry points:")
for label, indices in found_points.items():
    print(f"  {label}: k-point indices {indices}")
    for idx in indices:
        kpt = data.kpoints[idx]
        print(f"      k-point {idx}: ({kpt[0]:7.4f}, {kpt[1]:7.4f}, {kpt[2]:7.4f})")

