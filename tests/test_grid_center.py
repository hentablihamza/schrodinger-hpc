import os

import numpy as np
import pytest

import grid_center as gc
from conftest import EXAMPLE_DIR, PDB_1M17


def test_compute_grid_simple():
    coords = np.array(
        [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 4.0, 0.0], [0.0, 0.0, 6.0]]
    )
    box = gc.compute_grid(coords, outer_pad=10.0)
    assert box.center == pytest.approx((0.5, 1.0, 1.5))
    assert box.ligand_size == pytest.approx((2.0, 4.0, 6.0))
    assert box.outer_box == pytest.approx((12.0, 14.0, 16.0))
    assert box.n_atoms == 4


def test_empty_raises():
    with pytest.raises(ValueError):
        gc.compute_grid(np.empty((0, 3)))


@pytest.mark.skipif(not os.path.exists(PDB_1M17), reason="example PDB missing")
def test_erlotinib_centroid():
    lig = os.path.join(EXAMPLE_DIR, "ligands", "AQ4_A_999.pdb")
    if not os.path.exists(lig):
        pytest.skip("extracted ligand not present")
    box = gc.compute_grid_from_file(lig)
    assert box.n_atoms == 29
    # Centroid of erlotinib in 1M17 (reproducible to 2 dp).
    assert box.center[0] == pytest.approx(22.014, abs=0.05)
    assert box.center[1] == pytest.approx(0.253, abs=0.05)
    assert box.center[2] == pytest.approx(52.794, abs=0.05)
