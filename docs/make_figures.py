"""Generate the repo's documentation figures (workflow diagram + EGFR pocket).

Run:  python docs/make_figures.py
Outputs:  docs/workflow.png, docs/egfr_pocket.png
All numbers in the pocket figure come from the real 1M17 example via the
license-free tools, not from hard-coded values.
"""

import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "tools"))


def workflow_figure(out_path):
    steps = [
        ("SMILES\n(input.smi)", "#e3f2fd"),
        ("LigPrep\n(ligprep.sh)", "#bbdefb"),
        ("3D ligands\n(.maegz)", "#e3f2fd"),
        ("Glide grid\n(grid-gen.sh)", "#c8e6c9"),
        ("Glide SP dock\n(glide-array.sh)", "#a5d6a7"),
        ("Parse + rank\n(parse_glide.py)", "#fff9c4"),
    ]
    fig, ax = plt.subplots(figsize=(13, 3.2))
    ax.set_xlim(0, len(steps) * 2.2)
    ax.set_ylim(0, 3)
    ax.axis("off")

    box_w, box_h, y = 1.8, 1.3, 1.0
    centers = []
    for i, (label, color) in enumerate(steps):
        x = 0.3 + i * 2.2
        box = FancyBboxPatch(
            (x, y), box_w, box_h,
            boxstyle="round,pad=0.06,rounding_size=0.12",
            linewidth=1.5, edgecolor="#37474f", facecolor=color,
        )
        ax.add_patch(box)
        ax.text(x + box_w / 2, y + box_h / 2, label, ha="center", va="center",
                fontsize=10.5, fontweight="bold", color="#263238")
        centers.append(x + box_w)
        if i > 0:
            arrow = FancyArrowPatch(
                (centers[i - 1] - box_w + box_w, y + box_h / 2),
                (x, y + box_h / 2),
                arrowstyle="-|>", mutation_scale=18, linewidth=1.6, color="#455a64",
            )
            ax.add_patch(arrow)

    # Receptor prep branch feeding the grid.
    ax.text(0.3 + 3 * 2.2 + box_w / 2, y - 0.55,
            "receptor.pdb -> prepwizard.sh\n+ grid_center.py (pocket centroid)",
            ha="center", va="top", fontsize=8.5, color="#2e7d32", style="italic")
    ax.text(len(steps) * 1.1, 2.75,
            "Schrodinger on Hellbender HPC  (account/partition dkhf3-lab, module schrodinger/v2024-2)",
            ha="center", va="center", fontsize=11, fontweight="bold", color="#1a237e")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out_path}")


def pocket_figure(out_path):
    import extract_pocket as ep
    import grid_center as gc

    pdb = os.path.join(ROOT, "examples", "egfr", "1M17.pdb")
    if not os.path.exists(pdb):
        print("1M17.pdb missing; skipping pocket figure")
        return
    coords = ep.ligand_coords(pdb, ligand_res="AQ4")
    pocket, residues, _ = ep.find_pocket_residues(pdb, coords, cutoff=5.0)
    box = gc.compute_grid(coords)
    center = np.array(box.center)

    # Protein CA coords of pocket residues for a simple 3D-ish scatter.
    res_coords = []
    labels = []
    for res in residues:
        if "CA" in res:
            res_coords.append(res["CA"].get_coord())
            labels.append(f"{res.get_resname()}{res.id[1]}")
    res_coords = np.array(res_coords)

    fig = plt.figure(figsize=(8, 7))
    ax = fig.add_subplot(111, projection="3d")
    ax.scatter(coords[:, 0], coords[:, 1], coords[:, 2],
               c="#d32f2f", s=55, depthshade=True, label="Erlotinib (AQ4, 29 atoms)")
    ax.scatter(res_coords[:, 0], res_coords[:, 1], res_coords[:, 2],
               c="#1565c0", s=40, alpha=0.8,
               label=f"Pocket residues within 5 A ({len(pocket)})")
    ax.scatter(*center, c="#2e7d32", marker="*", s=380,
               edgecolor="k", label="Grid center (centroid)")
    for xyz, lab in zip(res_coords, labels):
        ax.text(xyz[0], xyz[1], xyz[2], lab, fontsize=6.5, color="#0d47a1")
    ax.set_title("EGFR ATP site (PDB 1M17) — grid center and 5 A pocket\n"
                 f"center = ({center[0]:.2f}, {center[1]:.2f}, {center[2]:.2f})",
                 fontsize=11)
    ax.set_xlabel("x (A)"); ax.set_ylabel("y (A)"); ax.set_zlabel("z (A)")
    ax.legend(loc="upper left", fontsize=8)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out_path}")


if __name__ == "__main__":
    workflow_figure(os.path.join(HERE, "workflow.png"))
    pocket_figure(os.path.join(HERE, "egfr_pocket.png"))
