#!/bin/bash
#SBATCH --error=mpi-n2-%j.err
#SBATCH --output=mpi-n2-%j.out
#SBATCH --partition=snowy
#SBATCH --time=01:00:00
#SBATCH --nodes=8
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=32G #SBATCH --job-name=mpi

module load gcc/8.3.0
module load openmpi/3.1.4

mpicxx -std=c++14 -fopenmp n2.cpp n2.h particle.h -o n2 -O3
mpirun -np 1 n2 body_1500.data > n2_1_1500.out
mpirun -np 2 n2 body_1500.data > n2_2_1500.out
mpirun -np 3 n2 body_1500.data > n2_3_1500.out
mpirun -np 4 n2 body_1500.data > n2_4_1500.out
mpirun -np 5 n2 body_1500.data > n2_5_1500.out
mpirun -np 6 n2 body_1500.data > n2_6_1500.out
mpirun -np 7 n2 body_1500.data > n2_7_1500.out
mpirun -np 8 n2 body_1500.data > n2_8_1500.out

# mpirun -np 1 n2 body_2000.data > n2_1_2000.out
# mpirun -np 7 n2 body_2000.data > n2_7_2000.out
# mpirun -np 8 n2 body_2000.data > n2_8_2000.out

# mpirun -np 1 n2 body_2500.data > n2_1_2500.out
# mpirun -np 7 n2 body_2500.data > n2_7_2500.out
# mpirun -np 8 n2 body_2500.data > n2_8_2500.out