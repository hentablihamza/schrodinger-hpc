#!/bin/bash
# LigPrep: 2D SMILES -> 3D, ionization/tautomer states (Epik), stereochemistry.
# Produces the .maegz ligand file that Glide docks. Real Hellbender header.
#
# Submit:   sbatch ligprep.sh
# Input:    input.smi (one SMILES + optional name per line)
# Output:   ligands.maegz

#SBATCH --job-name ligprep
#SBATCH --time 2-00:00:00
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 4
#SBATCH --mem-per-cpu=8gb
#SBATCH --account dkhf3-lab
#SBATCH --partition dkhf3-lab
#SBATCH --output ligprep-%j.out

module load schrodinger/v2024-2

# -ma 500   max 500 structures per input molecule (stereo/ionization/tautomer)
# -bff 14   OPLS_2005 force field
# -s 1      one low-energy ring conformation
# -epik     Epik for ionization/tautomer states at target pH
# -nt       neutralize before ionization
# -HOST localhost:4 matches --cpus-per-task above
ligprep -ma 500 -bff 14 -s 1 -epik -nt \
  -ismi input.smi -omae ligands.maegz \
  -NOJOBID -HOST localhost:4
