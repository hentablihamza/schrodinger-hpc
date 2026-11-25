"""Extract binding-pocket residues around a ligand from a PDB file.

Given a receptor PDB and a reference ligand (by residue name/chain, or an
external ligand PDB/SDF), find every protein residue with at least one atom
within a distance cutoff of the ligand and write them to a pocket PDB. The
residue list and the ligand centroid are what you feed to Glide grid
generation (grid center + box).

Schrodinger-free: Biopython for the receptor, optional RDKit for SDF ligands.
"""

from __future__ import annotations

import argparse
import os
from dataclasses import dataclass
from typing import Iterable

import numpy as np
from Bio.PDB import PDBParser, PDBIO, Select, NeighborSearch

from extract_ligands import AMINO_ACIDS


@dataclass
class PocketResidue:
    resname: str
    chain: str
    resseq: int
    icode: str
    min_dist: float

    @property
    def label(self) -> str:
        ic = self.icode.strip()
        return f"{self.resname}{self.resseq}{ic}:{self.chain}"


class _ResidueSetSelect(Select):
    def __init__(self, residues):
        self._keep = set(id(r) for r in residues)

    def accept_residue(self, residue):
        return 1 if id(residue) in self._keep else 0


def _ligand_coords_from_pdb(pdb_path, ligand_res, chain, structure_id="lig"):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(structure_id, pdb_path)
    coords = []
    for model in structure:
        for ch in model:
            if chain is not None and ch.id != chain:
                continue
            for residue in ch:
                rn = residue.get_resname().strip().upper()
                if ligand_res is not None and rn != ligand_res.upper():
                    continue
                hetflag = residue.id[0]
                if ligand_res is None and hetflag == " ":
                    continue
                for atom in residue:
                    coords.append(atom.get_coord())
        break
    return np.asarray(coords, dtype=float)


def _ligand_coords_from_sdf(sdf_path):
    from rdkit import Chem  # imported lazily so RDKit stays optional

    suppl = Chem.SDMolSupplier(sdf_path, removeHs=False)
    mol = next((m for m in suppl if m is not None), None)
    if mol is None:
        raise ValueError(f"No valid molecule in {sdf_path}")
    conf = mol.GetConformer()
    return np.asarray(
        [list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())],
        dtype=float,
    )


def ligand_coords(
    receptor_pdb: str,
    ligand_res: str | None = None,
    ligand_chain: str | None = None,
    ligand_file: str | None = None,
) -> np.ndarray:
    """Return the ligand's atom coordinates (N, 3)."""
    if ligand_file is not None:
        ext = os.path.splitext(ligand_file)[1].lower()
        if ext == ".sdf":
            return _ligand_coords_from_sdf(ligand_file)
        return _ligand_coords_from_pdb(ligand_file, ligand_res, ligand_chain)
    if ligand_res is None:
        raise ValueError("Provide either --ligand-res or --ligand-file")
    return _ligand_coords_from_pdb(receptor_pdb, ligand_res, ligand_chain)


def find_pocket_residues(
    receptor_pdb: str,
    lig_coords: np.ndarray,
    cutoff: float = 5.0,
    protein_only: bool = True,
    structure_id: str = "receptor",
):
    """Return (list[PocketResidue], Biopython residues, structure)."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(structure_id, receptor_pdb)
    model = next(iter(structure))

    atoms = list(model.get_atoms())
    ns = NeighborSearch(atoms)

    hits: dict[int, tuple] = {}
    for coord in lig_coords:
        for atom in ns.search(coord, cutoff, level="A"):
            residue = atom.get_parent()
            rn = residue.get_resname().strip().upper()
            if protein_only and rn not in AMINO_ACIDS:
                continue
            d = float(np.linalg.norm(atom.get_coord() - coord))
            key = id(residue)
            if key not in hits or d < hits[key][1]:
                hits[key] = (residue, d)

    residues = [r for r, _ in hits.values()]
    pocket = []
    for residue, dist in hits.values():
        chain = residue.get_parent().id
        hetflag, resseq, icode = residue.id
        pocket.append(
            PocketResidue(
                resname=residue.get_resname().strip().upper(),
                chain=chain,
                resseq=resseq,
                icode=icode,
                min_dist=round(dist, 2),
            )
        )
    pocket.sort(key=lambda p: (p.chain, p.resseq))
    return pocket, residues, structure


def write_pocket_pdb(structure, residues, out_path: str) -> str:
    io = PDBIO()
    io.set_structure(structure)
    io.save(out_path, _ResidueSetSelect(residues))
    return out_path


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description="Extract binding-pocket residues within a cutoff of a ligand."
    )
    ap.add_argument("receptor", help="receptor PDB file")
    ap.add_argument("--ligand-res", help="ligand residue name inside the receptor PDB")
    ap.add_argument("--ligand-chain", help="restrict ligand to this chain")
    ap.add_argument("--ligand-file", help="external ligand PDB or SDF")
    ap.add_argument("-c", "--cutoff", type=float, default=5.0, help="distance cutoff (A)")
    ap.add_argument("-o", "--out", default="pocket.pdb", help="output pocket PDB")
    ap.add_argument("--keep-het", action="store_true", help="keep non-protein residues too")
    args = ap.parse_args(argv)

    coords = ligand_coords(
        args.receptor, args.ligand_res, args.ligand_chain, args.ligand_file
    )
    pocket, residues, structure = find_pocket_residues(
        args.receptor, coords, cutoff=args.cutoff, protein_only=not args.keep_het
    )
    print(
        f"Ligand has {len(coords)} atoms; "
        f"{len(pocket)} residue(s) within {args.cutoff} A:"
    )
    for p in pocket:
        print(f"  {p.label:14s}  min dist {p.min_dist:.2f} A")
    write_pocket_pdb(structure, residues, args.out)
    print(f"Wrote pocket to {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
