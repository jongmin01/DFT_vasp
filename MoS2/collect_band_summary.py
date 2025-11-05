#!/usr/bin/env python3
# Robust EIGENVAL collector: scans all k-headers, then reads NBANDS lines if valid.
# Prints VBM at Γ and at K, plus (K-Γ). Comments are in English.

import re
from pathlib import Path

float_re = re.compile(r'[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?')
int_re   = re.compile(r'^[+-]?\d+$')

def is_k_header(line: str) -> bool:
    """k-header has at least 3 floats (kx, ky, kz) optionally a 4th (weight)."""
    parts = line.split()
    if len(parts) < 3: return False
    # First three tokens must be floats
    for t in parts[:3]:
        if not float_re.fullmatch(t): return False
    # ok
    return True

def is_band_line(line: str) -> bool:
    """Band line typically: integer_index energy [occupancy]."""
    parts = line.split()
    if len(parts) < 2: return False
    if not int_re.fullmatch(parts[0]): return False
    return bool(float_re.fullmatch(parts[1]))

def read_lines(p: Path):
    with p.open() as f:
        return f.readlines()

def parse_header_nk_nb(lines):
    """
    Find the 'nel nkpts nbands' triplet. In your files it looks like:
         24    360     96
    We'll search the first 120 lines for a line with exactly 3 integers.
    """
    for i in range(min(120, len(lines))):
        parts = lines[i].split()
        if len(parts) == 3 and all(int_re.fullmatch(x) for x in parts):
            try:
                nel, nk, nb = map(int, parts)
                return i, nel, nk, nb
            except:  # noqa
                pass
    # Fallback: accept lines with >=3 tokens; try last 3 as ints
    for i in range(min(120, len(lines))):
        parts = lines[i].split()
        if len(parts) >= 3 and all(int_re.fullmatch(x) for x in parts[-3:]):
            try:
                nel, nk, nb = map(int, parts[-3:])
                return i, nel, nk, nb
            except:  # noqa
                pass
    raise RuntimeError("Could not locate 'nel nkpts nbands' line")

def collect_k_headers(lines, start_idx, max_scan=100000):
    """Return list of indices where k-headers occur, scanning from start_idx+1 onward."""
    ks = []
    i = start_idx + 1
    end = min(len(lines), start_idx + 1 + max_scan)
    while i < end:
        # skip blank lines
        while i < end and lines[i].strip() == "":
            i += 1
        if i >= end: break
        if is_k_header(lines[i].strip()):
            ks.append(i)
            # jump ahead by a conservative step; actual advance will be done by reading bands
            i += 1
        else:
            i += 1
    return ks

def read_band_block(lines, k_idx, nbands):
    """
    Given index of k-header, read the subsequent band block (skip one optional blank).
    Return list of energies if exactly nbands lines are valid, else None.
    """
    off = k_idx + 1
    if off < len(lines) and lines[off].strip() == "":
        off += 1
    vals = []
    for j in range(nbands):
        idx = off + j
        if idx >= len(lines): return None
        ln = lines[idx]
        if not is_band_line(ln): return None
        parts = ln.split()
        try:
            e = float(parts[1])
        except Exception:
            return None
        vals.append(e)
    return vals

def vbm(energies):
    occ = [e for e in energies if e <= 0.0]
    return max(occ) if occ else None

def summarize_dir(d: Path):
    eig = d / "EIGENVAL"
    out = d / "OUTCAR"

    if not eig.exists() or eig.stat().st_size < 1000:
        return "EIGENVAL_TOO_SHORT"

    lines = read_lines(eig)
    hdr_idx, nel, nkpts, nbands = parse_header_nk_nb(lines)

    # Gather candidate k-headers then keep only those that have a valid band block
    k_candidates = collect_k_headers(lines, hdr_idx)
    k_blocks = []
    for k_idx in k_candidates:
        vals = read_band_block(lines, k_idx, nbands)
        if vals is not None:
            k_blocks.append((k_idx, vals))
        if len(k_blocks) >= nkpts:
            break

    if len(k_blocks) < nkpts:
        return f"K_BLOCKS_INCOMPLETE({len(k_blocks)}/{nkpts})"

    # Γ–K–M–Γ path: per-segment length ≈ NKPTS/3
    N = round(nkpts / 3)
    if N < 2:
        return f"BAD_SEGMENT_N(NKPTS={nkpts})"

    vG = vbm(k_blocks[0][1])
    vK = vbm(k_blocks[N-1][1])
    if vG is None or vK is None:
        return "NO_VBM_FOUND"
    dv = vK - vG
    return f"{vG:8.3f} {vK:12.3f} {dv:14.3f}"

def main():
    print("strain_dir, VBM(Gamma)[eV], VBM(K)[eV], (K-Gamma)[eV]  (K>Gamma → direct)")
    for d in sorted(Path(".").glob("strain_*/24_band_nosoc_HR")):
        res = summarize_dir(d)
        if " " in res and res.strip()[0] in "-0123456789":
            # numeric line
            print(f"{str(d.parent):<16s} {res}")
        else:
            print(f"{str(d.parent):<16s} {res}")

if __name__ == "__main__":
    main()
