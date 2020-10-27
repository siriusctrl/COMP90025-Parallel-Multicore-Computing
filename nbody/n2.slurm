#!/bin/bash
#SBATCH --error=mpi-n2-%j.err
#SBATCH --output=mpi-n2-%j.out
#SBATCH --partition=snowy
#SBATCH --time=01:00:00
#SBATCH --nodes=10
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=32G #SBATCH --job-name=mpi

module load gcc/8.3.0
module load openmpi/3.1.4

mpicxx -std=c++14 -fopenmp n2.cpp n2.h particle.h -o n2 -O3
mpirun -np 6 n2 body_1500.data