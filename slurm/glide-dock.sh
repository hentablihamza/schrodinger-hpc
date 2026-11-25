#!/bin/bash
# Glide SP docking of a prepared ligand library against a prebuilt grid.
# Serial job (Glide manages its own threading). Real University of Missouri
# "Hellbender" header, as used in production in the Kireev Lab.
#
# Submit:   sbatch glide-dock.sh
# Input:    an .in file (see egfr.in below) + the grid .zip + the .maegz ligands
# Output:   <name>_pv.maegz (poses) and, with glide_sort / -write_csv, a CSV.

#SBATCH --job-name egfr-glide
#SBATCH --time 3-00:00:00
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --mem-per-cpu=16gb
#SBATCH --account dkhf3-lab
#SBATCH --partition dkhf3-lab
#SBATCH --output glide-%j.out

module load schrodinger/v2024-2

# The .in file points to the grid and the prepared ligands, e.g.:
#   FORCEFIELD   OPLS_2005
#   GRIDFILE     glide-grid_egfr.zip
#   LIGANDFILE   egfr_ligands.maegz
#   PRECISION    SP
glide egfr.in -NOJOBID
