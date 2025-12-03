#!/usr/bin/env python3
"""
Fat-band plotting workflow using PyProcar with aggressive caching.
The first run parses PROCAR and stores ebs.pkl/structure.pkl. Afterwards
all plots are created straight from the binary cache, avoiding repeated
multi-gigabyte PROCAR reads.
"""

from __future__ import annotations

import gzip
import os
import shutil
from pathlib import Path
from typing import Iterable, List, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from pyprocar import io
from pyprocar.cfg import ConfigFactory, PlotType
from pyprocar.plotter import EBSPlot
from pyprocar.utils import data_utils

# Paths and caching helpers
WORKDIR = Path(__file__).resolve().parent
MPLCONFIG_DIR = WORKDIR / ".matplotlib"
PROCAR_PATH = WORKDIR / "PROCAR"
PROCAR_GZ_PATH = WORKDIR / "PROCAR.gz"
EBS_CACHE = WORKDIR / "ebs.pkl"
STRUCTURE_CACHE = WORKDIR / "structure.pkl"

# Energy window (eV, referenced to E_F)
ENERGY_MIN = -3.0
ENERGY_MAX = 3.0

# High-symmetry labels from Γ-M-K-Γ path
KPATH_LABELS = ["$\\Gamma$", "M", "K", "$\\Gamma$"]
KPATH_TICKS = [0, 43, 83, 120]

# Atom indices (0-based, extracted from POSCAR)
ATOMS_MO = list(range(0, 25))
ATOMS_W = list(range(25, 50))
ATOMS_S = list(range(50, 100))
ATOMS_SE = list(range(100, 150))

# Fatband styling knobs to improve readability
WEIGHT_MASK = 0.08  # drop contributions below 8%
WEIGHT_POWER = 0.6  # down-weight large markers (sqrt-ish)
LINEWIDTH_BASE = 0.28  # base linewidth before PyProcar's x5 boost
YTICKS = np.arange(ENERGY_MIN, ENERGY_MAX + 1, 1)


def ensure_matplotlib_cache() -> None:
    """Ensure matplotlib has a writable cache directory."""
    MPLCONFIG_DIR.mkdir(exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(MPLCONFIG_DIR))


def ensure_procar() -> None:
    """Ensure an uncompressed PROCAR is available for PyProcar parsing."""
    if PROCAR_PATH.exists():
        return
    if not PROCAR_GZ_PATH.exists():
        raise FileNotFoundError("Neither PROCAR nor PROCAR.gz is available")

    print("  PROCAR not found. Extracting PROCAR.gz (one-time operation)...")
    with gzip.open(PROCAR_GZ_PATH, "rb") as src, open(PROCAR_PATH, "wb") as dst:
        shutil.copyfileobj(src, dst)
    print("  Extraction complete.")


def load_cached_ebs():
    """Load cached ElectronicBandStructure data or parse it once from PROCAR."""
    ensure_matplotlib_cache()
    if EBS_CACHE.exists():
        print("  Loading cached PyProcar data (ebs.pkl)...")
        ebs = data_utils.load_pickle(EBS_CACHE)
        structure = (
            data_utils.load_pickle(STRUCTURE_CACHE) if STRUCTURE_CACHE.exists() else None
        )
        return ebs, structure

    ensure_procar()
    print("  Cache not found. Parsing PROCAR -> ebs.pkl (this may take a while)...")
    parser = io.Parser(code="vasp", dirpath=str(WORKDIR))
    ebs = parser.ebs
    structure = parser.structure
    data_utils.save_pickle(ebs, EBS_CACHE)
    data_utils.save_pickle(structure, STRUCTURE_CACHE)
    print("  Cached PyProcar data written to ebs.pkl/structure.pkl.")
    return ebs, structure


def build_config(title: str, cmap: str, *, show_colorbar: bool):
    """Create a PyProcar plotting config with shared styling."""
    config = ConfigFactory.create_config(PlotType.BAND_STRUCTURE)
    config.title = title
    config.cmap = cmap
    config.linewidth = [LINEWIDTH_BASE, LINEWIDTH_BASE]
    config.plot_color_bar = show_colorbar
    config.colorbar_title = "Projection weight"
    config.clim = [0.0, 1.0]
    config.figure_size = [7.0, 5.0]
    config.grid = True
    config.fermi_color = "#444444"
    return config


def finalize_plot(ebs_plot: EBSPlot, savefig: str) -> None:
    """Apply shared axis formatting, save, and close the figure."""
    ebs_plot.set_xticks(KPATH_TICKS, KPATH_LABELS)
    ebs_plot.set_yticks(interval=list(YTICKS))
    ebs_plot.set_ylim([ENERGY_MIN, ENERGY_MAX])
    ebs_plot.set_xlabel()
    ebs_plot.set_ylabel()
    ebs_plot.draw_fermi(fermi_level=0.0)
    ebs_plot.grid()
    ebs_plot.set_title()
    ebs_plot.save(savefig)
    plt.close(ebs_plot.fig)


def plot_total_bands(ebs, spins: Sequence[int]) -> str:
    """Plot the plain band structure."""
    print("\n1. Plotting total band structure...")
    config = build_config("Total Band Structure", cmap="viridis", show_colorbar=False)
    ebs_plot = EBSPlot(ebs, ebs.kpath, spins=spins, config=config)
    ebs_plot.plot_bands()
    finalize_plot(ebs_plot, "fatband_total.png")
    print("   Saved: fatband_total.png")
    return "fatband_total.png"


def _prepare_weights(raw_weights: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """
    Normalize and soften weights for more legible fatbands.
    Returns (width_weights, color_weights).
    """
    weights = np.clip(raw_weights, 0.0, None)
    max_w = float(weights.max()) + 1e-12
    scaled = np.power(weights / max_w, WEIGHT_POWER)
    return scaled, scaled


def plot_parametric(
    ebs,
    spins: Sequence[int],
    atoms: Iterable[int],
    *,
    orbitals: Iterable[int] | None,
    title: str,
    cmap: str,
    output_name: str,
) -> str:
    """Generic helper for element/orbital fat-band projections."""
    weights_raw = ebs.ebs_sum(
        atoms=list(atoms),
        orbitals=None if orbitals is None else list(orbitals),
        spins=list(spins),
    )
    width_weights, color_weights = _prepare_weights(weights_raw)
    config = build_config(title, cmap=cmap, show_colorbar=True)
    ebs_plot = EBSPlot(ebs, ebs.kpath, spins=spins, config=config)
    ebs_plot.plot_parameteric(
        spins=list(spins),
        width_weights=width_weights,
        color_weights=color_weights,
        width_mask=WEIGHT_MASK,
        color_mask=None,
        elimit=[ENERGY_MIN, ENERGY_MAX],
        labels=[title],
    )
    ebs_plot.set_colorbar_title(title=f"{title} weight")
    finalize_plot(ebs_plot, output_name)
    print(f"   Saved: {output_name}")
    return output_name


def plot_fatband_by_element(ebs, spins: Sequence[int]) -> List[str]:
    """Produce total and element-projected fat bands."""
    print("=" * 60)
    print("Generating fat-band plots using cached PyProcar data")
    print("=" * 60)
    generated = [plot_total_bands(ebs, spins)]

    element_jobs = [
        ("Band Structure - Mo Contribution", ATOMS_MO, "Reds", "fatband_Mo.png"),
        ("Band Structure - W Contribution", ATOMS_W, "Blues", "fatband_W.png"),
        ("Band Structure - S Contribution", ATOMS_S, "Greens", "fatband_S.png"),
        ("Band Structure - Se Contribution", ATOMS_SE, "Purples", "fatband_Se.png"),
    ]
    for title, atoms, cmap, filename in element_jobs:
        print(f"\nPlotting {title}...")
        generated.append(
            plot_parametric(
                ebs,
                spins,
                atoms,
                orbitals=None,
                title=title,
                cmap=cmap,
                output_name=filename,
            )
        )
    return generated


def plot_orbital_fatband(ebs, spins: Sequence[int]) -> List[str]:
    """Produce orbital-resolved fat bands for selected species."""
    print("\n" + "=" * 60)
    print("Generating orbital-resolved fat-band plots")
    print("=" * 60)

    outputs = []
    orbitals_d = [4, 5, 6, 7, 8]  # dxy, dyz, dz2, dxz, dx2
    orbitals_p = [1, 2, 3]  # py, pz, px

    orbital_jobs = [
        ("Band Structure - Mo d orbitals", ATOMS_MO, orbitals_d, "Reds", "fatband_Mo_d.png"),
        ("Band Structure - W d orbitals", ATOMS_W, orbitals_d, "Blues", "fatband_W_d.png"),
        ("Band Structure - S p orbitals", ATOMS_S, orbitals_p, "Greens", "fatband_S_p.png"),
        ("Band Structure - Se p orbitals", ATOMS_SE, orbitals_p, "Purples", "fatband_Se_p.png"),
    ]

    for title, atoms, orbitals, cmap, filename in orbital_jobs:
        print(f"\nPlotting {title}...")
        outputs.append(
            plot_parametric(
                ebs,
                spins,
                atoms,
                orbitals=orbitals,
                title=title,
                cmap=cmap,
                output_name=filename,
            )
        )
    return outputs


def main():
    ebs, _ = load_cached_ebs()
    spins = list(range(ebs.nspins))
    generated_files = []
    generated_files.extend(plot_fatband_by_element(ebs, spins))
    generated_files.extend(plot_orbital_fatband(ebs, spins))

    print("\n" + "=" * 60)
    print("Fat-band generation complete. Files created:")
    for fname in generated_files:
        print(f"  - {fname}")
    print("=" * 60)


if __name__ == "__main__":
    main()
