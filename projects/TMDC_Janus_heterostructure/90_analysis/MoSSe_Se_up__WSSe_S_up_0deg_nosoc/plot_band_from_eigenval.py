#!/usr/bin/env python3
"""
Plot band structure directly from an EIGENVAL file.

Reads EIGENVAL for k-points/energies, DOSCAR for E_F, detects high-symmetry
breaks along the path, and writes a simple band plot shifted by the Fermi level.
"""

import argparse
from typing import List, Tuple

import matplotlib.pyplot as plt
import numpy as np

# Reuse readers already present in this folder
from analyze_bandgap_eigenval import read_doscar_fermi, read_eigenval


def kpath_distance(kpoints: np.ndarray) -> np.ndarray:
    """Return cumulative distance along the k-path."""
    distances = [0.0]
    for i in range(1, len(kpoints)):
        dk = np.linalg.norm(kpoints[i] - kpoints[i - 1])
        distances.append(distances[-1] + dk)
    return np.array(distances)


def label_for_kpoint(kpt: np.ndarray, tol: float = 1e-3, fallback: str = "") -> str:
    """Assign a simple label to common hexagonal high-symmetry points."""
    if np.allclose(kpt, [0, 0, 0], atol=tol):
        return "Gamma"
    if (np.allclose(kpt, [1 / 3, 1 / 3, 0], atol=tol) or
            np.allclose(kpt, [-1 / 3, -1 / 3, 0], atol=tol) or
            np.allclose(kpt, [2 / 3, 2 / 3, 0], atol=tol)):
        return "K"
    if (np.allclose(kpt, [0.5, 0, 0], atol=tol) or
            np.allclose(kpt, [0, 0.5, 0], atol=tol) or
            np.allclose(kpt, [0.5, 0.5, 0], atol=tol)):
        return "M"
    return fallback


def detect_high_symmetry_points(
    kpoints: np.ndarray,
    angle_threshold: float = 0.995,
    tol: float = 1e-6
) -> Tuple[List[int], List[str]]:
    """
    Detect path break points using direction changes or repeated k-points.

    Returns indices for vertical guide lines and their labels.
    """
    candidates: List[int] = [0]

    for i in range(1, len(kpoints) - 1):
        step_prev = kpoints[i] - kpoints[i - 1]
        step_next = kpoints[i + 1] - kpoints[i]
        prev_norm = np.linalg.norm(step_prev)
        next_norm = np.linalg.norm(step_next)

        # Zero-length step means an intentional repeat between path segments
        if prev_norm < tol or next_norm < tol:
            candidates.append(i)
            continue

        cos_angle = float(np.dot(step_prev, step_next) / (prev_norm * next_norm))
        if cos_angle < angle_threshold:
            candidates.append(i)

    candidates.append(len(kpoints) - 1)

    seen = set()
    indices: List[int] = []
    labels: List[str] = []
    for idx in candidates:
        if idx in seen:
            continue
        seen.add(idx)
        indices.append(idx)
        labels.append(label_for_kpoint(kpoints[idx], fallback=f"k{idx+1}"))

    return indices, labels


def plot_bandstructure(
    kpoints: np.ndarray,
    energies: np.ndarray,
    efermi: float,
    output: str,
    ymin: float,
    ymax: float
) -> None:
    """Plot the band structure and save to disk."""
    kdist = kpath_distance(kpoints)
    energies_shifted = energies - efermi
    high_sym_indices, high_sym_labels = detect_high_symmetry_points(kpoints)

    fig, ax = plt.subplots(figsize=(8, 6))

    for iband in range(energies.shape[1]):
        ax.plot(kdist, energies_shifted[:, iband], color="black", linewidth=0.8, alpha=0.8)

    ax.axhline(0.0, color="red", linestyle="--", linewidth=1.0, alpha=0.8, label="E_F")

    for idx in high_sym_indices:
        ax.axvline(kdist[idx], color="gray", linestyle="-", linewidth=0.6, alpha=0.5)

    ax.set_xticks([kdist[i] for i in high_sym_indices])
    ax.set_xticklabels(high_sym_labels)

    ax.set_ylabel("Energy - E_F (eV)")
    ax.set_xlabel("k-path")
    ax.set_xlim(kdist[0], kdist[-1])
    ax.set_ylim(ymin, ymax)
    ax.legend(loc="upper right", frameon=False)
    ax.grid(axis="y", linestyle="--", alpha=0.3)

    fig.tight_layout()
    fig.savefig(output, dpi=300)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot band structure from EIGENVAL.")
    parser.add_argument("--eigenval", default="EIGENVAL", help="Path to EIGENVAL file.")
    parser.add_argument("--doscar", default="DOSCAR", help="Path to DOSCAR file for Fermi level.")
    parser.add_argument("--output", default="bandstructure_eigenval.png", help="Output image name.")
    parser.add_argument(
        "--window",
        type=float,
        default=3.0,
        help="Energy window (eV) plotted around E_F when not using full range or custom bounds.",
    )
    parser.add_argument("--emin", type=float, default=None, help="Lower bound relative to E_F (eV).")
    parser.add_argument("--emax", type=float, default=None, help="Upper bound relative to E_F (eV).")
    parser.add_argument(
        "--full-range",
        action="store_true",
        help="Plot the full energy span found in EIGENVAL instead of a window around E_F.",
    )
    args = parser.parse_args()

    kpoints, energies, nkpts, nbands = read_eigenval(args.eigenval)
    efermi = read_doscar_fermi(args.doscar)

    energies_shifted = energies - efermi
    if args.full_range:
        ymin = energies_shifted.min() - 0.25
        ymax = energies_shifted.max() + 0.25
    elif args.emin is not None or args.emax is not None:
        ymin = args.emin if args.emin is not None else energies_shifted.min() - 0.25
        ymax = args.emax if args.emax is not None else energies_shifted.max() + 0.25
    else:
        ymin = -abs(args.window)
        ymax = abs(args.window)

    # Ensure sane limits
    if ymax <= ymin:
        span = max(1.0, abs(ymin))
        ymax = ymin + span

    plot_bandstructure(kpoints, energies, efermi, args.output, ymin, ymax)
    print(f"Plotted {nbands} bands over {nkpts} k-points to {args.output}")


if __name__ == "__main__":
    main()
