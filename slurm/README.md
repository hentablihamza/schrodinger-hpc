# SLURM job scripts

Ready-to-adapt SLURM batch scripts for running Schrödinger on the University of
Missouri **Hellbender** HPC cluster. The working values baked in are the real
ones used in production in the Kireev Lab:

- **Account & partition:** `dkhf3-lab` (both the SLURM account and the partition)
- **Module:** `module load schrodinger/v2024-2`

> Kireev is the **lab / PI**, not a partition. The partition and account are
> both `dkhf3-lab`.

## The scripts

| Script | What it does | Key SLURM settings |
|---|---|---|
| `ligprep.sh` | 2D SMILES → 3D ligands with ionization/tautomer states (Epik) | `--cpus-per-task 4`, `--time 2-00:00:00` |
| `prepwizard.sh` | Protein Preparation Wizard (add H, fill side chains, protonation) | `--cpus-per-task 2`, `--time 06:00:00` |
| `grid-gen.sh` | Build a Glide receptor grid centered on the binding site | serial, `--time 02:00:00` |
| `glide-dock.sh` | Glide SP docking of one ligand file against a grid | serial, `--mem-per-cpu=16gb`, `--time 3-00:00:00` |
| `glide-array.sh` | **Job array** for large / ultra-large virtual screening | `--array=0-199%40` |
| `longleaf-glide-gpu.sh` | **Historical**, second cluster: UNC Longleaf GPU variant | `--partition volta-gpu`, `--gres=gpu:1` |

## Workflow order

```
input.smi ──ligprep.sh──▶ ligands.maegz
receptor.pdb ──prepwizard.sh──▶ receptor_prepared.maegz
         + grid_center.py (centroid) ──grid-gen.sh──▶ glide-grid_<target>.zip
ligands.maegz + grid ──glide-dock.sh / glide-array.sh──▶ *_pv.maegz + CSV
         ──parse_glide.py──▶ ranked CSV
```

## Submitting and monitoring

```bash
sbatch ligprep.sh                 # submit
squeue -u "$USER"                 # your queued/running jobs
sacct -j <jobid> --format=JobID,State,Elapsed,MaxRSS   # accounting

# Summarize a big array job's progress (license-free, offline-testable):
sacct -j <jobid> --format=JobID,State,Elapsed -P | python ../tools/squeue_summary.py -
```

## Large-scale screening with `glide-array.sh`

The array script is the workhorse for screening large libraries. Split a
prepared library into `chunks/chunk_000.maegz … chunk_NNN.maegz` first, then:

```bash
sbatch glide-array.sh     # --array=0-199%40  => 200 chunks, 40 at a time
```

Each array task reads `chunks/chunk_${SLURM_ARRAY_TASK_ID}.maegz`, writes its
own `.in`, and docks only that chunk. Failed indices are trivially
re-submittable (`sbatch --array=3,57 glide-array.sh`).

## Adapting to your own cluster

Change these to your site's values:

- `--account <your-account>` and `--partition <your-partition>`
- `module load schrodinger/<your-version>`
- walltime / memory to match your scheduler's limits

Everything else (the Schrödinger commands and the `.in` file format) is
cluster-independent.

## A note on CI

Schrödinger is **licensed commercial software** and is never run in CI. These
scripts are syntax-checked (`bash -n`) and linted with **ShellCheck** on every
push; the Schrödinger binaries themselves run only on the licensed cluster.
