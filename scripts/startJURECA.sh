#!/bin/bash -l
#SBATCH --nodes=16
#SBATCH --ntasks=512
#SBATCH --ntasks-per-node=32
#SBATCH --time=12:00:00
#SBATCH --output=mpi-out.%j
#SBATCH --error=mpi_err.%j
#SBATCH --partition=batch 

srun ../bin/ifos3d  ./in_and_out/ifos3d_toy.inp | tee ./in_and_out/ifos3d_toy.out

