#!/bin/bash
#SBATCH --job-name=[job name]
#SBATCH --nodes=200
#SBATCH --tasks-per-node=8
#SBATCH --cpus-per-task=8
#SBATCH --account=[account name]
#SBATCH --partition=standard
#SBATCH --qos=standard
#SBATCH --time=00:00:00

module load namd/2.14
export OMP_NUM_THREADS=8

srun --unbuffered --cpu-bind=cores namd2 +setcpuaffinity +ppn 7 protein_input.namd > simulation_out.log
