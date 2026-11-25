"""Reproduce the EGFR (PDB 1M17) worked example end to end, license-free.

Runs the real helper tools on 1M17 and prints a short report. This is what
produced the committed outputs (ligands/, pocket.pdb, grid_center.json) and the
docs figure. No Schrodinger needed.

Run:  python examples/egfr/run_example.py
"""

import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(ROOT, "tools"))

import extract_ligands as el
import extract_pocket as ep
import grid_center as gc
import parse_glide as pg

PDB = os.path.join(HERE, "1M17.pdb")


def main():
    # 1. Extract ligands (erlotinib = AQ4; waters are skipped).
    ligands = el.extract_all(PDB, os.path.join(HERE, "ligands"))
    print(f"[1] Extracted {len(ligands)} ligand(s):")
    for lig, path in ligands:
        print(f"    {lig.label}  ({lig.n_atoms} atoms) -> {os.path.basename(path)}")

    # 2. Binding-pocket residues within 5 A of erlotinib.
    coords = ep.ligand_coords(PDB, ligand_res="AQ4")
    pocket, residues, structure = ep.find_pocket_residues(PDB, coords, cutoff=5.0)
    ep.write_pocket_pdb(structure, residues, os.path.join(HERE, "pocket.pdb"))
    print(f"[2] {len(pocket)} pocket residues within 5 A:")
    print("    " + ", ".join(p.label.split(":")[0] for p in pocket))

    # 3. Grid center + box from the ligand centroid.
    box = gc.compute_grid(coords)
    with open(os.path.join(HERE, "grid_center.json"), "w") as fh:
        json.dump(box.as_dict(), fh, indent=2)
    cx, cy, cz = box.center
    print(f"[3] Grid center: {cx:.3f}, {cy:.3f}, {cz:.3f}")

    # 4. Parse the (SYNTHETIC) example Glide CSV to show the ranking step.
    ranked = pg.parse_glide_csv(os.path.join(HERE, "glide_results_example.csv"))
    ranked.to_csv(os.path.join(HERE, "glide_results_ranked.csv"), index=False)
    print(f"[4] Parsed example Glide CSV -> {len(ranked)} ranked ligands "
          f"(top hit: {ranked.iloc[0]['title']} @ {ranked.iloc[0]['docking_score']})")
    print("    NOTE: the Glide CSV is SYNTHETIC/illustrative; Glide itself needs a license.")


if __name__ == "__main__":
    main()
