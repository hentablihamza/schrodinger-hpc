#!/bin/bash
# Protein Preparation Wizard: add hydrogens, assign bond orders, fill in missing
# side chains, set protonation/tautomer states (Epik/PROPKA), and restrained
# minimization. Run this on the raw PDB before building a grid.
#
# Submit:   sbatch prepwizard.sh
# Input:    receptor.pdb (or .mae)
# Output:   receptor_prepared.maegz

#SBATCH --job-name prepwizard
#SBATCH --time 06:00:00
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 2
#SBATCH --mem-per-cpu=8gb
#SBATCH --account dkhf3-lab
#SBATCH --partition dkhf3-lab
#SBATCH --output prepwizard-%j.out

module load schrodinger/v2024-2

# -fillsidechains  rebuild missing side-chain atoms
# -disulfides      detect disulfide bonds
# -epik_pH 7.0     protonation states at physiological pH
# -propka_pH 7.0   PROPKA pKa prediction
# -minimize_adj_h  optimize hydrogen-bond network
"${SCHRODINGER}/utilities/prepwizard" \
  -fillsidechains -disulfides \
  -epik_pH 7.0 -propka_pH 7.0 -minimize_adj_h \
  receptor.pdb receptor_prepared.maegz \
  -NOJOBID -HOST localhost:2
