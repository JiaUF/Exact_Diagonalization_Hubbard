#!/bin/bash

#SBATCH --job-name=AA_10_hf
#SBATCH --output=status.%J.out
#SBATCH --error=status.%J.err

#SBATCH --time=02:00:00
#SBATCH --partition=iric,simes,owners
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=40G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jiacjustc@gmail.com

module load math arpack
module load devel icc
icc -o main main.cc arpack.cc hamiltonian.cc -larpack -lopenblas

srun -n 1 main
