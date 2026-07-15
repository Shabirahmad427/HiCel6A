#!/bin/bash
#SBATCH -n 16                         # Number of CPU cores (adjust as needed)
#SBATCH --job-name=pcasuite
#SBATCH --nodes=1
#SBATCH --partition=par48-i          # or par48-x depending on your user profile
#SBATCH -o par48out.txt
#SBATCH -e par48err.txt
#SBATCH --mail-type=ALL
#SBATCH --time=72:00:00

cd $SLURM_SUBMIT_DIR

# === Environment setup ===
module load hwloc
module load gnu12
module load openmpi4   # harmless even if pcazip doesn’t use MPI

export OMP_NUM_THREADS=$SLURM_NTASKS

# === Commands ===
echo "Running PCAzip on $SLURM_NTASKS cores"
pcazip -i protein.traj.nc -p protein.avg.pdb -o protein.traj.pcz -q 100

# Dump first three principal component projections as PDB animations
pczdump -i protein.traj.pcz --anim 1 --pdb > protein.traj.proj1.pdb
pczdump -i protein.traj.pcz --anim 2 --pdb > protein.traj.proj2.pdb
pczdump -i protein.traj.pcz --anim 3 --pdb > protein.traj.proj3.pdb

echo "Done."

