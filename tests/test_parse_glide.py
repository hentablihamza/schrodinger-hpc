import os

import pandas as pd
import pytest

import parse_glide as pg
from conftest import EXAMPLE_DIR

CSV = os.path.join(EXAMPLE_DIR, "glide_results_example.csv")


@pytest.mark.skipif(not os.path.exists(CSV), reason="example CSV missing")
def test_best_per_ligand_ranks():
    df = pg.parse_glide_csv(CSV, best_per_ligand=True)
    # 10 unique ligands; duplicate erlotinib/gefitinib poses collapsed.
    assert len(df) == 10
    assert list(df["rank"]) == list(range(1, 11))
    # Best score is the most negative, and the top hit is erlotinib.
    assert df.iloc[0]["title"] == "erlotinib"
    assert df.iloc[0]["docking_score"] == pytest.approx(-9.82)
    # Monotonic non-decreasing (less negative) down the ranking.
    scores = df["docking_score"].tolist()
    assert scores == sorted(scores)


@pytest.mark.skipif(not os.path.exists(CSV), reason="example CSV missing")
def test_all_poses_kept():
    df = pg.parse_glide_csv(CSV, best_per_ligand=False)
    assert len(df) == 12  # two extra poses in the example file


@pytest.mark.skipif(not os.path.exists(CSV), reason="example CSV missing")
def test_top_n():
    df = pg.parse_glide_csv(CSV, top=3)
    assert len(df) == 3
    assert df.iloc[-1]["rank"] == 3


def test_missing_score_column(tmp_path):
    p = tmp_path / "bad.csv"
    pd.DataFrame({"title": ["x"], "foo": [1]}).to_csv(p, index=False)
    with pytest.raises(ValueError):
        pg.parse_glide_csv(str(p))
