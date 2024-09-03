#!/bin/bash

#SBATCH --qos=cjia1-b
#SBATCH --job-name=1DHubbard
#SBATCH --output=status.%J.out
#SBATCH --error=status.%J.err

#SBATCH --time=0:59:59
##SBATCH --partition=bigmem
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jiacjustc@gmail.com

module load mkl
module load intel
icc -o main main.cc arpack.cc hamiltonian.cc -mkl
srun -n 1 main

