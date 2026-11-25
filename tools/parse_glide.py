"""Parse Glide docking results into a tidy, ranked table.

Glide writes poses to ``*_pv.maegz`` and (with ``-write_csv`` / via
``glide_sort``) a CSV with per-pose properties. The ``.maegz`` pose file needs
Schrodinger to read, so this tool works from the **CSV export**, which is plain
text and license-free. It finds the GlideScore / docking-score column
(``r_i_docking_score`` or ``r_i_glide_gscore``), keeps the best pose per ligand
(or all poses), ranks by score, and writes a clean CSV.

Schrodinger-free: pandas only.
"""

from __future__ import annotations

import argparse

import pandas as pd

# Candidate column names Glide uses for the main docking score, best first.
SCORE_COLUMNS = [
    "r_i_docking_score",
    "r_i_glide_gscore",
    "docking_score",
    "GlideScore",
    "glide_gscore",
    "r_i_glide_emodel",  # last resort (emodel, not gscore)
]

# Candidate columns identifying the ligand/title.
TITLE_COLUMNS = ["title", "s_m_title", "Title", "NAME", "name", "ligand", "Ligand"]


def _pick(columns, candidates):
    lower = {c.lower(): c for c in columns}
    for cand in candidates:
        if cand in columns:
            return cand
        if cand.lower() in lower:
            return lower[cand.lower()]
    return None


def parse_glide_csv(
    path: str,
    best_per_ligand: bool = True,
    top: int | None = None,
) -> pd.DataFrame:
    """Read a Glide CSV export and return a ranked DataFrame.

    Columns: ``title``, ``docking_score``, ``rank`` (plus the original score
    column name preserved). Lower (more negative) GlideScore is better.
    """
    df = pd.read_csv(path)
    score_col = _pick(df.columns, SCORE_COLUMNS)
    if score_col is None:
        raise ValueError(
            f"No recognizable Glide score column in {path}. "
            f"Looked for: {', '.join(SCORE_COLUMNS)}. Found: {list(df.columns)}"
        )
    title_col = _pick(df.columns, TITLE_COLUMNS)

    out = pd.DataFrame()
    out["title"] = df[title_col] if title_col else [f"pose_{i+1}" for i in range(len(df))]
    out["docking_score"] = pd.to_numeric(df[score_col], errors="coerce")
    out = out.dropna(subset=["docking_score"])

    if best_per_ligand:
        out = (
            out.sort_values("docking_score")
            .groupby("title", as_index=False, sort=False)
            .first()
        )

    out = out.sort_values("docking_score").reset_index(drop=True)
    out.insert(0, "rank", range(1, len(out) + 1))
    if top is not None:
        out = out.head(top).reset_index(drop=True)
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Parse and rank Glide CSV docking results.")
    ap.add_argument("csv", help="Glide CSV export")
    ap.add_argument("-o", "--out", help="write ranked CSV here")
    ap.add_argument("--top", type=int, help="keep only the top N ligands")
    ap.add_argument(
        "--all-poses", action="store_true",
        help="keep every pose instead of the best per ligand",
    )
    args = ap.parse_args(argv)

    ranked = parse_glide_csv(
        args.csv, best_per_ligand=not args.all_poses, top=args.top
    )
    if args.out:
        ranked.to_csv(args.out, index=False)
        print(f"Wrote {len(ranked)} ranked rows to {args.out}")
    with pd.option_context("display.max_rows", 20):
        print(ranked.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
