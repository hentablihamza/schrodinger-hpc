#!/bin/bash
# Large-scale / ultra-large virtual screening as a SLURM JOB ARRAY.
#
# Split a big prepared ligand library into N equal chunks ahead of time
# (chunks/chunk_000.maegz ... chunk_NNN.maegz, plus one .in per chunk, or a
# template .in that references $CHUNK). Each array task docks exactly one chunk,
# so hundreds of chunks dock in parallel across the cluster, bounded by %K.
#
#   --array=0-199%40  ->  200 chunks, at most 40 running at once.
#
# This is the workhorse pattern for screening millions of compounds: one tiny
# serial Glide job per chunk, scheduled independently, trivially restartable
# (just resubmit the failed array indices).
#
# Submit:   sbatch glide-array.sh
# Monitor:  sacct -j <jobid> --format=JobID,State,Elapsed -P | python ../tools/squeue_summary.py -

#SBATCH --job-name egfr-vs-array
#SBATCH --time 1-00:00:00
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --mem-per-cpu=16gb
#SBATCH --account dkhf3-lab
#SBATCH --partition dkhf3-lab
#SBATCH --array=0-199%40
#SBATCH --output logs/vs-%A_%a.out

module load schrodinger/v2024-2

mkdir -p results logs

# Zero-pad the array index to match chunk file names (chunk_000.maegz, ...).
CHUNK=$(printf "%03d" "${SLURM_ARRAY_TASK_ID}")
LIGANDS="chunks/chunk_${CHUNK}.maegz"
INFILE="results/chunk_${CHUNK}.in"

if [[ ! -f "${LIGANDS}" ]]; then
  echo "No ligand chunk ${LIGANDS}; nothing to do for task ${SLURM_ARRAY_TASK_ID}."
  exit 0
fi

# Build this task's Glide input from the shared grid + this chunk of ligands.
cat > "${INFILE}" <<INEOF
FORCEFIELD	OPLS_2005
GRIDFILE	glide-grid_egfr.zip
LIGANDFILE	${LIGANDS}
PRECISION	SP
INEOF

cd results || exit 1
glide "chunk_${CHUNK}.in" -NOJOBID
