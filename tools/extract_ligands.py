"""Extract HETATM ligands from a PDB file.

A crystal structure often contains many non-protein groups: the ligand you
care about, but also crystallographic waters, ions, and buffer/cryoprotectant
molecules. Before you can build a docking grid or prepare a receptor you need
to pull out the *real* ligand(s) and drop the rest.

This module is Schrodinger-free: it only needs Biopython, so it runs anywhere
(including CI and a login node without a license).
"""

from __future__ import annotations

import argparse
import os
from dataclasses import dataclass, field
from typing import Iterable

from Bio.PDB import PDBParser, PDBIO, Select

# Residue names that are almost never the ligand of interest. We skip these so
# the extracted set is the genuine small-molecule ligand(s).
WATERS = {"HOH", "WAT", "DOD", "H2O"}

IONS = {
    "NA", "K", "MG", "CA", "ZN", "MN", "FE", "FE2", "CU", "CU1", "NI", "CO",
    "CD", "HG", "CL", "BR", "IOD", "F", "SO4", "PO4", "PI", "CAC", "NH4",
}

# Common crystallization additives / cryoprotectants / buffers.
BUFFERS = {
    "GOL", "EDO", "PEG", "PG4", "PGE", "1PE", "2PE", "DMS", "MPD", "ACT",
    "ACY", "FMT", "TRS", "EPE", "MES", "BME", "DTT", "IMD", "CIT", "TLA",
    "MLI", "SCN", "AZI", "NO3", "BCT", "OXL", "FLC", "P6G", "7PE", "12P",
    "15P", "PE4", "SIN", "SRT", "GBL", "BOG", "LDA", "C8E", "OLC",
}

# Standard amino acids and nucleotides are polymer residues, not ligands.
AMINO_ACIDS = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    "MSE", "SEC", "PYL", "HID", "HIE", "HIP", "CYX", "ASH", "GLH", "LYN",
}

NUCLEOTIDES = {"A", "C", "G", "U", "T", "DA", "DC", "DG", "DT", "DU", "N"}

DEFAULT_SKIP = WATERS | IONS | BUFFERS | AMINO_ACIDS | NUCLEOTIDES


@dataclass
class Ligand:
    """A single extracted ligand (one HETATM residue)."""

    resname: str
    chain: str
    resseq: int
    icode: str
    n_atoms: int
    residue: object = field(default=None, repr=False)

    @property
    def label(self) -> str:
        ic = self.icode.strip()
        return f"{self.resname}_{self.chain}_{self.resseq}{ic}"


class _LigandSelect(Select):
    def __init__(self, residue):
        self._res = residue

    def accept_residue(self, residue):
        return 1 if residue is self._res else 0


def _is_candidate(residue, skip: Iterable[str], min_atoms: int) -> bool:
    hetflag, _, _ = residue.id
    # hetflag is " " for standard polymer residues, "W" for water, "H_XXX" for
    # heteroatoms. We only want genuine HETATM groups.
    if hetflag == " ":
        return False
    resname = residue.get_resname().strip().upper()
    if resname in skip:
        return False
    n_atoms = sum(1 for _ in residue.get_atoms())
    if n_atoms < min_atoms:
        return False
    return True


def find_ligands(
    pdb_path: str,
    skip: Iterable[str] | None = None,
    min_atoms: int = 6,
    structure_id: str = "structure",
) -> list[Ligand]:
    """Return every HETATM ligand in a PDB, skipping waters/ions/buffers.

    ``min_atoms`` filters out tiny groups (single atoms, small ions that slip
    past the name list). Raise it to be stricter about what counts as a ligand.
    """

    skip_set = {s.strip().upper() for s in (skip if skip is not None else DEFAULT_SKIP)}
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(structure_id, pdb_path)

    ligands: list[Ligand] = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if not _is_candidate(residue, skip_set, min_atoms):
                    continue
                hetflag, resseq, icode = residue.id
                ligands.append(
                    Ligand(
                        resname=residue.get_resname().strip().upper(),
                        chain=chain.id,
                        resseq=resseq,
                        icode=icode,
                        n_atoms=sum(1 for _ in residue.get_atoms()),
                        residue=residue,
                    )
                )
        break  # first model only
    return ligands


def write_ligand_pdb(ligand: Ligand, out_path: str) -> str:
    """Write a single ligand to its own PDB file."""
    # Rebuild a parent structure containing only this residue so PDBIO works.
    structure = ligand.residue.get_parent().get_parent().get_parent()
    io = PDBIO()
    io.set_structure(structure)
    io.save(out_path, _LigandSelect(ligand.residue))
    return out_path


def extract_all(
    pdb_path: str,
    out_dir: str,
    skip: Iterable[str] | None = None,
    min_atoms: int = 6,
) -> list[tuple[Ligand, str]]:
    """Extract every ligand to ``out_dir`` as ``<label>.pdb``. Returns pairs."""
    os.makedirs(out_dir, exist_ok=True)
    ligands = find_ligands(pdb_path, skip=skip, min_atoms=min_atoms)
    written: list[tuple[Ligand, str]] = []
    for lig in ligands:
        out_path = os.path.join(out_dir, f"{lig.label}.pdb")
        write_ligand_pdb(lig, out_path)
        written.append((lig, out_path))
    return written


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Extract HETATM ligands from a PDB file.")
    ap.add_argument("pdb", help="input PDB file")
    ap.add_argument("-o", "--out-dir", default="ligands", help="output directory")
    ap.add_argument(
        "--min-atoms", type=int, default=6,
        help="ignore HETATM groups with fewer atoms (default: 6)",
    )
    ap.add_argument(
        "--keep-all", action="store_true",
        help="do not skip waters/ions/buffers",
    )
    ap.add_argument(
        "--list", action="store_true",
        help="only list ligands, do not write files",
    )
    args = ap.parse_args(argv)

    skip = set() if args.keep_all else None
    if args.list:
        ligands = find_ligands(args.pdb, skip=skip, min_atoms=args.min_atoms)
        print(f"Found {len(ligands)} ligand(s) in {args.pdb}:")
        for lig in ligands:
            print(f"  {lig.label:20s}  {lig.n_atoms:3d} atoms  (resname {lig.resname})")
        return 0

    written = extract_all(args.pdb, args.out_dir, skip=skip, min_atoms=args.min_atoms)
    print(f"Extracted {len(written)} ligand(s) to {args.out_dir}/:")
    for lig, path in written:
        print(f"  {os.path.basename(path):24s}  {lig.n_atoms:3d} atoms")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
