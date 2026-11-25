#!/bin/bash
# HISTORICAL / SECOND CLUSTER: UNC "Longleaf" GPU variant.
#
# This is kept as a reference example of running Glide on a *different* cluster
# (UNC Longleaf) with a GPU partition and an older Schrodinger release. It is
# NOT the current Hellbender workflow -- the active scripts above target
# Hellbender (account/partition dkhf3-lab, module schrodinger/v2024-2).
# Shown to document multi-cluster experience; adapt paths to your own site.

#SBATCH --job-name glide-gpu
#SBATCH --time 3-00:00:00
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --mem-per-cpu=16gb
#SBATCH --partition volta-gpu
#SBATCH --qos gpu_access
#SBATCH --gres=gpu:1
#SBATCH --output glide-gpu-%j.out

module load schrodinger/2021-4
module load cuda/9.2

# On Longleaf the Glide binary lives under the cluster apps tree.
/nas/longleaf/apps/schrodinger/2021-4/glide egfr.in -NOJOBID
