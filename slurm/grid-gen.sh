#!/bin/bash
# Generate a Glide receptor grid from a prepared receptor, centered on the
# binding site. Get the GRID_CENTER coordinates from the co-crystal ligand with
# the license-free helper:
#
#   python ../tools/grid_center.py ligand_erlotinib.pdb
#
# then paste the centroid into the GRID_CENTER line of the .in file below.
#
# Submit:   sbatch grid-gen.sh
# Input:    receptor_prepared.maegz  (+ the grid .in written below)
# Output:   glide-grid_egfr.zip

#SBATCH --job-name grid-gen
#SBATCH --time 02:00:00
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --mem-per-cpu=8gb
#SBATCH --account dkhf3-lab
#SBATCH --partition dkhf3-lab
#SBATCH --output grid-%j.out

module load schrodinger/v2024-2

# GRID_CENTER is the ligand centroid from tools/grid_center.py.
# INNERBOX / OUTERBOX are the ligand-diameter and search boxes (Angstrom).
cat > grid.in <<INEOF
GRIDFILE	glide-grid_egfr.zip
RECEP_FILE	receptor_prepared.maegz
GRID_CENTER	22.014, 0.253, 52.794
INNERBOX	10, 10, 10
OUTERBOX	30, 30, 30
INEOF

glide grid.in -NOJOBID
