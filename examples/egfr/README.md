# Worked example: EGFR (PDB 1M17, erlotinib)

A fully license-free walk-through of the *preparation and analysis* steps that
surround a Glide docking run, using the EGFR kinase domain co-crystallized with
erlotinib (PDB **1M17**, ligand residue **AQ4**). Every number here is produced
by running the helper tools on the real structure — reproduce it with:

```bash
python examples/egfr/run_example.py
```

## What the tools find (real outputs)

**1. Ligand extraction** (`tools/extract_ligands.py`)
- 1 ligand extracted: **AQ4 (erlotinib), 29 atoms**, chain A → `ligands/AQ4_A_999.pdb`
- 20 crystallographic waters correctly skipped

**2. Binding-pocket residues** (`tools/extract_pocket.py`, 5 Å cutoff)
- **20 residues** around erlotinib → `pocket.pdb`:
  LEU694, GLY695, VAL702, ALA719, **LYS721** (catalytic Lys), GLU738, MET742,
  LEU764, ILE765, **THR766** (gatekeeper), GLN767, LEU768, **MET769** (hinge),
  PRO770, PHE771, GLY772, CYS773, LEU820, THR830, **ASP831** (DFG Asp)
- These are the canonical EGFR ATP-site residues — a good sanity check that the
  pocket extraction is correct.

**3. Grid center** (`tools/grid_center.py`)
- Centroid of erlotinib: **(22.014, 0.253, 52.794)** → `grid_center.json`
- This is the `GRID_CENTER` value already pasted into `slurm/grid-gen.sh`.

**4. Result parsing** (`tools/parse_glide.py`)
- `glide_results_example.csv` is **SYNTHETIC / illustrative** (Glide itself needs
  a license, so no real docking scores are fabricated here). It shows the parser
  collapsing duplicate poses to the best per ligand and ranking them →
  `glide_results_ranked.csv`.

## How these feed the SLURM scripts

```
1M17.pdb ──extract_pocket.py──▶ pocket residues (sanity check)
         ──extract_ligands.py──▶ AQ4_A_999.pdb
         ──grid_center.py──────▶ GRID_CENTER (22.014, 0.253, 52.794)
                                      │
receptor_prepared.maegz + GRID_CENTER ──grid-gen.sh──▶ glide-grid_egfr.zip
your_library.smi ──ligprep.sh──▶ ligands.maegz
ligands.maegz + grid ──glide-array.sh──▶ poses + CSV ──parse_glide.py──▶ ranked hits
```

## Files

| File | Source | Notes |
|---|---|---|
| `1M17.pdb` | RCSB PDB | Real crystal structure (attribution below) |
| `ligands/AQ4_A_999.pdb` | generated | Extracted erlotinib |
| `pocket.pdb` | generated | 20 pocket residues (5 Å) |
| `grid_center.json` | generated | Centroid + box |
| `glide_results_example.csv` | **synthetic** | Illustrates the parser only |
| `glide_results_ranked.csv` | generated | Parser output |
| `sacct_array_sample.txt` | synthetic | Sample SLURM output for `squeue_summary.py` |

## Data attribution

Structure **1M17** from the RCSB Protein Data Bank (Stamos, Sliwkowski &
Eigenbrot, *J. Biol. Chem.* 2002, 277:46265). PDB data are freely available from
[rcsb.org](https://www.rcsb.org/).
