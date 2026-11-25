"""Compute a docking grid center and box size from a ligand.

Glide (and every other docking engine) needs a grid box: a center point and
an edge length. The most reliable way to place it on a known binding site is
to take the co-crystal ligand, use its centroid as the center, and size the
box from the ligand's extent plus padding.

Schrodinger-free: RDKit for SDF/MOL ligands, a tiny PDB reader otherwise.
"""

from __future__ import annotations

import argparse
import os
from dataclasses import dataclass

import numpy as np


@dataclass
class GridBox:
    center: tuple[float, float, float]
    inner_box: tuple[float, float, float]   # Glide inner box (ligand diameter)
    outer_box: tuple[float, float, float]   # Glide outer box (search space)
    n_atoms: int
    ligand_size: tuple[float, float, float]  # bounding-box extents

    def as_dict(self) -> dict:
        return {
            "center_x": round(self.center[0], 3),
            "center_y": round(self.center[1], 3),
            "center_z": round(self.center[2], 3),
            "inner_x": round(self.inner_box[0], 1),
            "inner_y": round(self.inner_box[1], 1),
            "inner_z": round(self.inner_box[2], 1),
            "outer_x": round(self.outer_box[0], 1),
            "outer_y": round(self.outer_box[1], 1),
            "outer_z": round(self.outer_box[2], 1),
            "n_atoms": self.n_atoms,
        }


def _coords_from_pdb(path: str) -> np.ndarray:
    coords = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(("ATOM", "HETATM")):
                coords.append(
                    (float(line[30:38]), float(line[38:46]), float(line[46:54]))
                )
    return np.asarray(coords, dtype=float)


def _coords_from_rdkit(path: str) -> np.ndarray:
    from rdkit import Chem

    ext = os.path.splitext(path)[1].lower()
    if ext == ".sdf":
        suppl = Chem.SDMolSupplier(path, removeHs=False)
        mol = next((m for m in suppl if m is not None), None)
    else:
        mol = Chem.MolFromMolFile(path, removeHs=False)
    if mol is None:
        raise ValueError(f"Could not read a molecule from {path}")
    conf = mol.GetConformer()
    return np.asarray(
        [list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())],
        dtype=float,
    )


def ligand_coords(path: str) -> np.ndarray:
    ext = os.path.splitext(path)[1].lower()
    if ext in (".sdf", ".mol"):
        return _coords_from_rdkit(path)
    return _coords_from_pdb(path)


def compute_grid(
    coords: np.ndarray,
    inner_pad: float = 0.0,
    outer_pad: float = 10.0,
) -> GridBox:
    """Centroid + Glide-style inner/outer boxes.

    Inner box ~ ligand bounding box (where ligand centers may sit); outer box
    adds padding (the full search space). Defaults follow common Glide practice
    (outer = ligand extent + ~10 A, giving roughly a 20-30 A cubic box).
    """
    if coords.size == 0:
        raise ValueError("No coordinates provided")
    center = coords.mean(axis=0)
    mins = coords.min(axis=0)
    maxs = coords.max(axis=0)
    extent = maxs - mins
    inner = extent + inner_pad
    outer = extent + outer_pad
    return GridBox(
        center=tuple(center),
        inner_box=tuple(inner),
        outer_box=tuple(outer),
        n_atoms=int(coords.shape[0]),
        ligand_size=tuple(extent),
    )


def compute_grid_from_file(path: str, **kw) -> GridBox:
    return compute_grid(ligand_coords(path), **kw)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description="Compute grid center and box size from a ligand (PDB/SDF/MOL)."
    )
    ap.add_argument("ligand", help="ligand file (PDB, SDF, or MOL)")
    ap.add_argument("--inner-pad", type=float, default=0.0, help="inner box padding (A)")
    ap.add_argument("--outer-pad", type=float, default=10.0, help="outer box padding (A)")
    ap.add_argument("--json", action="store_true", help="emit JSON")
    args = ap.parse_args(argv)

    box = compute_grid_from_file(
        args.ligand, inner_pad=args.inner_pad, outer_pad=args.outer_pad
    )
    if args.json:
        import json

        print(json.dumps(box.as_dict(), indent=2))
        return 0

    cx, cy, cz = box.center
    print(f"Ligand atoms:   {box.n_atoms}")
    print(f"Grid center:    {cx:8.3f} {cy:8.3f} {cz:8.3f}")
    print(
        "Ligand extent:  "
        f"{box.ligand_size[0]:.1f} x {box.ligand_size[1]:.1f} x {box.ligand_size[2]:.1f} A"
    )
    print(
        "Inner box:      "
        f"{box.inner_box[0]:.1f} x {box.inner_box[1]:.1f} x {box.inner_box[2]:.1f} A"
    )
    print(
        "Outer box:      "
        f"{box.outer_box[0]:.1f} x {box.outer_box[1]:.1f} x {box.outer_box[2]:.1f} A"
    )
    print()
    print("Feed the center into Glide grid generation (GRID_CENTER) in grid-gen.sh.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
