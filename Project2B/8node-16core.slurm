#!/bin/bash

#SBATCH --error=mpi-%j.err
#SBATCH --output=mpi-%j.out
#SBATCH --partition=snowy
#SBATCH --time=0:40:00
#SBATCH --nodes=8-8
#SBATCH --cpus-per-task=16
#SBATCH --ntasks-per-node=1
#SBATCH --mem=32G #SBATCH --job-name=mpi 
# You need to load a compiler before openmpi.
module load gcc/8.3.0
module load openmpi/3.1.4

mpicxx -std=c++14 -fopenmp -o p2 p2.cpp -O3
mpirun -np 8 p3 < mseq-big13-example.dat