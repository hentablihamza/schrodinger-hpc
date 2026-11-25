import os

import numpy as np
import pytest

import extract_pocket as ep
from conftest import PDB_1M17


pytestmark = pytest.mark.skipif(
    not os.path.exists(PDB_1M17), reason="1M17.pdb example not present"
)


def test_ligand_coords():
    coords = ep.ligand_coords(PDB_1M17, ligand_res="AQ4")
    assert coords.shape == (29, 3)


def test_pocket_residues_contain_hinge():
    coords = ep.ligand_coords(PDB_1M17, ligand_res="AQ4")
    pocket, residues, _ = ep.find_pocket_residues(PDB_1M17, coords, cutoff=5.0)
    labels = {p.label.split(":")[0] for p in pocket}
    # Known EGFR ATP-site residues around erlotinib.
    assert "MET769" in labels   # hinge
    assert "THR766" in labels   # gatekeeper
    assert "LYS721" in labels   # catalytic lysine
    assert 15 <= len(pocket) <= 30
    # All reported residues are genuinely within the cutoff.
    assert all(p.min_dist <= 5.0 for p in pocket)


def test_cutoff_monotonic():
    coords = ep.ligand_coords(PDB_1M17, ligand_res="AQ4")
    small, _, _ = ep.find_pocket_residues(PDB_1M17, coords, cutoff=4.0)
    big, _, _ = ep.find_pocket_residues(PDB_1M17, coords, cutoff=6.0)
    assert len(big) >= len(small)


def test_write_pocket(tmp_path):
    coords = ep.ligand_coords(PDB_1M17, ligand_res="AQ4")
    pocket, residues, structure = ep.find_pocket_residues(PDB_1M17, coords)
    out = tmp_path / "pocket.pdb"
    ep.write_pocket_pdb(structure, residues, str(out))
    assert out.exists()
