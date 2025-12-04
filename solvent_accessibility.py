#!/usr/bin/env python3
"""
优化后的溶剂可及性分析工具
通过只构建一次KD-tree来优化性能，提升约10%
"""

# 使用示例:python solvent_accessibility.py --wet-pdb SUMO1_water.pdb --dry-pdb SUMO1.pdb --method peratom --threshold 3.5 --margin 2.0 --fraction-threshold 0.3 --min-hits 2 --small-residue-size 5 --verbose


#!/usr/bin/env python3
import argparse
import csv
import os
import sys
import time
from collections import namedtuple  

import numpy as np
from Bio.PDB import PDBParser

# optional libs
try:
    from scipy.spatial import KDTree
    HAS_SCIPY = True
except Exception:
    HAS_SCIPY = False

ResidueInfo = namedtuple("ResidueInfo", ["chain", "resnum", "resname", "coord"])

# Default small-residue set
DEFAULT_SMALL_RES = {"GLY", "ALA", "SER", "THR", "CYS", "PRO"}

# ------------------------
# PDB parsing
# ------------------------
def parse_structure_get_residue_coords(pdb_file):
    parser = PDBParser(QUIET=True)
    struct = parser.get_structure("prot", pdb_file)

    residues = []
    water_coords = []
    water_names = {"HOH", "WAT", "SOL", "H2O", "TIP3", "TIP3P", "T3P", "W"}

    for model in struct:
        for chain in model:
            for residue in chain:
                rn = residue.get_resname().upper().strip()
                het = residue.id[0]

                # water detection
                if rn in water_names:
                    for atom in residue:
                        elem = getattr(atom, "element", "").upper()
                        aname = atom.get_name().strip().upper()
                        if aname in ("O", "OW", "OH", "OW1", "O1") or elem == "O":
                            water_coords.append(atom.coord)
                    continue

                # skip hetero (ligands etc.)
                if het.strip():
                    continue

                coords = [atom.coord for atom in residue
                          if getattr(atom, "element", "").upper() != "H"
                          and not atom.get_name().strip().startswith("H")]
                if not coords:
                    continue
                centroid = np.mean(np.array(coords, dtype=float), axis=0)
                residues.append(ResidueInfo(chain.id, residue.id[1], rn, centroid))

    return residues, np.array(water_coords, dtype=float), struct

# ------------------------
# centroid min distances (vectorized chunked)
# ------------------------
def compute_min_distances_chunked(res_coords, water_coords, chunk=5000):
    n_res = len(res_coords)
    if water_coords.size == 0:
        return np.full(n_res, np.inf)
    min_d2 = np.full(n_res, np.inf)
    n_w = len(water_coords)
    for start in range(0, n_w, chunk):
        end = min(start + chunk, n_w)
        diff = res_coords[:, None, :] - water_coords[None, start:end, :]
        d2 = np.sum(diff * diff, axis=2)
        min_d2 = np.minimum(min_d2, np.min(d2, axis=1))
    return np.sqrt(min_d2)

# ------------------------
# count waters within R
# ------------------------
def count_waters_within_R(res_coords, water_coords=None, R=5.0, chunk=5000, tree=None):
    if len(res_coords) == 0:
        return np.array([], dtype=int)
    
    if tree is not None:
        return np.array([len(tree.query_ball_point(rc, R)) for rc in res_coords], dtype=int)
    
    if water_coords is None or water_coords.size == 0:
        return np.zeros(len(res_coords), dtype=int)
    
    if HAS_SCIPY:
        tree = KDTree(water_coords)
        return np.array([len(tree.query_ball_point(rc, R)) for rc in res_coords], dtype=int)
    
    R2 = R * R
    counts = np.zeros(len(res_coords), dtype=int)
    n_w = len(water_coords)
    for start in range(0, n_w, chunk):
        end = min(start + chunk, n_w)
        diff = res_coords[:, None, :] - water_coords[None, start:end, :]
        d2 = np.sum(diff * diff, axis=2)
        counts += np.sum(d2 <= R2, axis=1)
    return counts

# ------------------------
# per-atom distances
# ------------------------
def collect_peratom_dists(struct, residues, water_coords=None, tree=None):
    dists_map = {}
    min_d_map = {}
    n_atoms_map = {}

    if tree is None:
        if HAS_SCIPY and water_coords is not None and water_coords.size != 0:
            tree = KDTree(water_coords)
        else:
            tree = None

    for r in residues:
        try:
            residue = struct[0][r.chain][(" ", r.resnum, " ")]
        except Exception:
            key = (r.chain, str(r.resnum))
            dists_map[key] = np.array([np.inf])
            min_d_map[key] = np.inf
            n_atoms_map[key] = 0
            continue

        atom_coords = []
        for atom in residue:
            if getattr(atom, "element", "").upper() == "H":
                continue
            if atom.get_name().strip().startswith("H"):
                continue
            atom_coords.append(atom.coord)
        if len(atom_coords) == 0:
            key = (r.chain, str(r.resnum))
            dists_map[key] = np.array([np.inf])
            min_d_map[key] = np.inf
            n_atoms_map[key] = 0
            continue

        atom_coords = np.array(atom_coords, dtype=float)
        if tree is not None:
            dists, _ = tree.query(atom_coords, k=1)
            dists = np.array(dists, dtype=float)
        else:
            if water_coords.size == 0:
                dists = np.full(len(atom_coords), np.inf)
            else:
                diff = atom_coords[:, None, :] - water_coords[None, :, :]
                d2 = np.sum(diff * diff, axis=2)
                dists = np.sqrt(np.min(d2, axis=1))
        key = (r.chain, str(r.resnum))
        dists_map[key] = dists
        min_d_map[key] = float(np.min(dists))
        n_atoms_map[key] = len(dists)
    return dists_map, n_atoms_map, min_d_map

# ------------------------
# CSV writers
# ------------------------
def write_csv(path, header, rows):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        for r in rows:
            w.writerow(r)

# ------------------------
# per-atom rules
# ------------------------
def peratom_flags_from_dists(residues, dists_map, fraction_threshold, min_hits, small_residue_size, small_res_set, threshold):
    flags = []
    min_d_list = []
    for r in residues:
        key = (r.chain, str(r.resnum))
        dists = dists_map.get(key, np.array([np.inf], dtype=float))
        n_atoms = len(dists)
        min_d = float(np.min(dists)) if n_atoms > 0 else np.inf
        min_d_list.append(min_d)

        is_small = (r.resname.upper() in small_res_set) or (n_atoms <= small_residue_size)
        n_hits = int((dists <= threshold).sum()) if n_atoms > 0 else 0
        fraction = float(n_hits) / float(n_atoms) if n_atoms > 0 else 0.0

        if is_small:
            accessible = (n_hits >= min_hits)
        else:
            accessible = (fraction >= fraction_threshold) and (n_hits >= min_hits)
        flags.append("Yes" if accessible else "No")
    return flags, np.array(min_d_list, dtype=float)

# -------------------------------------------------------------
# Wrap the original main into a function (do NOT change logic)
# -------------------------------------------------------------
def run_custom_accessibility(pdb, method, threshold, margin, R, chunk, nproc,
                             fraction_threshold, min_hits, small_residue_size,
                             verbose):

    prefix = os.path.basename(pdb).rsplit(".", 1)[0]
    start = time.time()

    residues, water_coords, struct = parse_structure_get_residue_coords(pdb)
    if verbose:
        print(f"Parsed {len(residues)} residues and {water_coords.shape[0]} water oxygens")

    res_coords = np.vstack([r.coord for r in residues]) if len(residues) > 0 else np.zeros((0, 3))
    min_d_centroid = compute_min_distances_chunked(res_coords, water_coords, chunk=chunk)

    water_tree = None
    if HAS_SCIPY and water_coords.size != 0:
        water_tree = KDTree(water_coords)
    
    counts = count_waters_within_R(res_coords, water_coords=water_coords, R=R, chunk=chunk, tree=water_tree)

    out_centroid = f"{prefix}_centroid.csv"
    centroid_rows = []
    for r, d, c in zip(residues, min_d_centroid, counts):
        centroid_rows.append([r.chain, r.resnum, r.resname, f"{d:.3f}", int(c), 
                              "Yes" if d <= threshold else "No"])
    write_csv(out_centroid,
              ["chain", "resnum", "resname", "minDist_A", "nWaterWithinR",
               f"accessible_{threshold}A"], centroid_rows)

    if verbose:
        acc = sum(1 for r in centroid_rows if r[-1] == "Yes")
        print(f"Wrote {out_centroid} Accessible={acc}/{len(residues)}")

    # per-atom
    dists_map, n_atoms_map, min_d_map = collect_peratom_dists(struct, residues,
                                                              water_coords=water_coords,
                                                              tree=water_tree)

    if method == "peratom":
        flags_peratom, min_d_final = peratom_flags_from_dists(
            residues, dists_map, fraction_threshold, min_hits,
            small_residue_size, DEFAULT_SMALL_RES, threshold
        )

        for i, (r, d_cent) in enumerate(zip(residues, min_d_centroid)):
            if d_cent > (threshold + margin):
                flags_peratom[i] = "No"
                min_d_final[i] = min_d_centroid[i]

    else:
        min_d_final = np.array([
            min_d_map.get((r.chain, str(r.resnum)), min_d_centroid[i])
            if (r.chain, str(r.resnum)) in min_d_map else min_d_centroid[i]
            for i, r in enumerate(residues)
        ])
        flags_peratom = ["Yes" if d <= threshold else "No" for d in min_d_final]

    out_peratom = f"{prefix}_peratom.csv"
    rows_peratom = []
    for idx, (r, d_key) in enumerate(zip(residues, min_d_final)):
        key = (r.chain, str(r.resnum))
        n_atoms = n_atoms_map.get(key, 0)
        nwater = int(counts[idx]) if len(counts) > 0 else 0
        flag = "Yes" if flags_peratom[idx] == "Yes" else "No"
        rows_peratom.append([r.chain, r.resnum, r.resname,
                             f"{float(d_key):.3f}", nwater, flag])
    write_csv(out_peratom,
              ["chain", "resnum", "resname", "minDist_A", "nWaterWithinR",
               f"accessible_{threshold}A"],
              rows_peratom)

    if verbose:
        acc = sum(1 for r in rows_peratom if r[-1] == "Yes")
        print(f"Wrote {out_peratom} Accessible={acc}/{len(residues)}")

    if verbose:
        print(f"Total runtime: {time.time() - start:.2f}s")

    # Return for comparison
    return residues, flags_peratom

# =====================================================================
# ------------------------ FreeSASA part -------------------------------
# =====================================================================
import freesasa

ACCESS_THRESHOLD = 10.0
WATER_NAMES = {"HOH", "WAT", "SOL", "H2O", "TIP3", "TIP3P", "T3P", "W"}

def compute_residue_sasa(pdb_file):
    structure = freesasa.Structure(pdb_file)
    result = freesasa.calc(structure)

    res_dict = result.residueAreas()
    out = []

    for chain, chain_dict in res_dict.items():
        for resnum, area_obj in chain_dict.items():
            resname = area_obj.residueType

            if resname.upper() in WATER_NAMES:
                continue

            sasa = area_obj.total
            accessible = "Yes" if sasa >= ACCESS_THRESHOLD else "No"
            out.append([chain, str(resnum), resname, sasa, accessible])
    return out

# =====================================================================
# ----------------------------- NEW MAIN ------------------------------
# =====================================================================

def main():
    parser = argparse.ArgumentParser(description="Combined solvent-accessibility + FreeSASA comparison")
    parser.add_argument("--wet-pdb", required=True, help="solvated pdb (for custom method)")
    parser.add_argument("--dry-pdb", required=True, help="unsolvated pdb (for FreeSASA)")
    parser.add_argument("--method", choices=("centroid","peratom"), default="peratom")
    parser.add_argument("--threshold", type=float, default=3.5)
    parser.add_argument("--margin", type=float, default=2.0)
    parser.add_argument("--R", type=float, default=5.0)
    parser.add_argument("--chunk", type=int, default=5000)
    parser.add_argument("--nproc", type=int, default=1)
    parser.add_argument("--fraction-threshold", type=float, default=0.20)
    parser.add_argument("--min-hits", type=int, default=1)
    parser.add_argument("--small-residue-size", type=int, default=5)
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    # ---------------------------------------------------------
    # Step 1. Run your original custom accessibility
    # ---------------------------------------------------------
    residues, custom_flags = run_custom_accessibility(
        args.wet_pdb,
        args.method,
        args.threshold,
        args.margin,
        args.R,
        args.chunk,
        args.nproc,
        args.fraction_threshold,
        args.min_hits,
        args.small_residue_size,
        args.verbose
    )

    # ---------------------------------------------------------
    # Step 2. Run FreeSASA on dry pdb
    # ---------------------------------------------------------
    sasa_results = compute_residue_sasa(args.dry_pdb)
    out_freesasa = args.dry_pdb.replace(".pdb", "_freesasa.csv")
    write_csv(out_freesasa,
              ["chain","resnum","resname","SASA","Accessible"],
              sasa_results)

    # ---------------------------------------------------------
    # Step 3. Compare (chain, resnum)
    # ---------------------------------------------------------
    def normalize_chain(c):
        c = c.strip()
        return c if c else "A"   # 如果是空链，统一当作 A（按你的实际情况可改）

    sasa_map = {(normalize_chain(c), str(r)): acc
            for c, r, rn, sasa, acc in sasa_results}
    custom_map = {(normalize_chain(r.chain), str(r.resnum)): custom_flags[i]
              for i, r in enumerate(residues)}
    
    comparison = []
    match_count = 0
    total = 0

    for r in residues:
        key = (normalize_chain(r.chain), str(r.resnum))
        c_acc = custom_map.get(key, "No")
        s_acc = sasa_map.get(key, "No")

        match = "Match" if c_acc == s_acc else "Mismatch"
        if match == "Match":
            match_count += 1
        total += 1

        comparison.append([
            r.chain, r.resnum, r.resname,
            c_acc, s_acc, match
        ])

    ratio = match_count / total if total > 0 else 0.0
    comparison.append(["","","","","",""])
    comparison.append(["Match_Ratio", f"{ratio:.4f}"])

    out_cmp = "comparison.csv"
    write_csv(out_cmp,
              ["chain","resnum","resname","Custom","FreeSASA","Match"],
              comparison)

    print("\n=== Comparison Completed ===")
    print(f"Generated: {out_freesasa}")
    print(f"Generated: comparison.csv")
    print(f"Match Ratio = {ratio:.4f}\n")


if __name__ == "__main__":
    main()
