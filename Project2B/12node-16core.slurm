#!/bin/bash

#SBATCH --error=mpi-12n16c-%j.err
#SBATCH --output=mpi-12n16c-%j.out
#SBATCH --partition=snowy
#SBATCH --time=0:10:00
#SBATCH --nodes=12
#SBATCH --cpus-per-task=16
#SBATCH --ntasks-per-node=1
#SBATCH --mem=32G #SBATCH --job-name=mpi 
# You need to load a compiler before openmpi.
module load gcc/8.3.0
module load openmpi/3.1.4

mpicxx -std=c++14 -fopenmp -o p3 p3.cpp -O3
mpirun -np 12 p3 < mseq-big13-example.dat
