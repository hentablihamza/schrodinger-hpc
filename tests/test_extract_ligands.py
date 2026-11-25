import os

import pytest

import extract_ligands as el
from conftest import PDB_1M17


pytestmark = pytest.mark.skipif(
    not os.path.exists(PDB_1M17), reason="1M17.pdb example not present"
)


def test_finds_erlotinib_only():
    ligands = el.find_ligands(PDB_1M17)
    # 1M17 holds erlotinib (AQ4) plus waters; waters must be skipped.
    assert len(ligands) == 1
    lig = ligands[0]
    assert lig.resname == "AQ4"
    assert lig.n_atoms == 29
    assert lig.chain == "A"


def test_skip_set_filters_waters():
    # With no skipping and a tiny min_atoms, waters (HOH) reappear.
    ligands = el.find_ligands(PDB_1M17, skip=set(), min_atoms=1)
    names = {l.resname for l in ligands}
    assert "HOH" in names
    assert "AQ4" in names


def test_write_ligand(tmp_path):
    written = el.extract_all(PDB_1M17, str(tmp_path))
    assert len(written) == 1
    lig, path = written[0]
    assert os.path.exists(path)
    with open(path) as fh:
        text = fh.read()
    assert "AQ4" in text
    # Only the ligand, no protein ATOM records of other residues.
    assert "HETATM" in text


def test_label():
    ligands = el.find_ligands(PDB_1M17)
    assert ligands[0].label == "AQ4_A_999"
