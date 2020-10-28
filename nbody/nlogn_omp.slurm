#!/bin/bash
#SBATCH --error=mpi-n2-%j.err
#SBATCH --output=mpi-n2-%j.out
#SBATCH --partition=snowy
#SBATCH --time=01:00:00
#SBATCH --nodes=4
#SBATCH --cpus-per-task=8
#SBATCH --ntasks-per-node=1
#SBATCH --mem=32G #SBATCH --job-name=mpi

module load gcc/8.3.0
module load openmpi/3.1.4

mpicxx -std=c++14 -fopenmp -o nlogn nlogn.h cell.h nlogn.cpp particle.h -O3
mpirun -np 4 nlogn body_1500.data 1 > nlogn_omp1_1500.out
mpirun -np 4 nlogn body_1500.data 2 > nlogn_omp2_1500.out
mpirun -np 4 nlogn body_1500.data 3 > nlogn_omp3_1500.out
mpirun -np 4 nlogn body_1500.data 4 > nlogn_omp4_1500.out
mpirun -np 4 nlogn body_1500.data 5 > nlogn_omp5_1500.out
mpirun -np 4 nlogn body_1500.data 6 > nlogn_omp6_1500.out
mpirun -np 4 nlogn body_1500.data 7 > nlogn_omp7_1500.out
mpirun -np 4 nlogn body_1500.data 8 > nlogn_omp8_1500.out

mpirun -np 4 nlogn body_1000.data 1 > nlogn_omp1_1000.out
mpirun -np 4 nlogn body_1000.data 2 > nlogn_omp2_1000.out
mpirun -np 4 nlogn body_1000.data 3 > nlogn_omp3_1000.out
mpirun -np 4 nlogn body_1000.data 4 > nlogn_omp4_1000.out
mpirun -np 4 nlogn body_1000.data 5 > nlogn_omp5_1000.out
mpirun -np 4 nlogn body_1000.data 6 > nlogn_omp6_1000.out
mpirun -np 4 nlogn body_1000.data 7 > nlogn_omp7_1000.out
mpirun -np 4 nlogn body_1000.data 8 > nlogn_omp8_1000.out

mpirun -np 4 nlogn body_2500.data 1 > nlogn_omp1_2500.out
mpirun -np 4 nlogn body_2500.data 2 > nlogn_omp2_2500.out
mpirun -np 4 nlogn body_2500.data 3 > nlogn_omp3_2500.out
mpirun -np 4 nlogn body_2500.data 4 > nlogn_omp4_2500.out
mpirun -np 4 nlogn body_2500.data 5 > nlogn_omp5_2500.out
mpirun -np 4 nlogn body_2500.data 6 > nlogn_omp6_2500.out
mpirun -np 4 nlogn body_2500.data 7 > nlogn_omp7_2500.out
mpirun -np 4 nlogn body_2500.data 8 > nlogn_omp8_2500.out

# mpirun -np 1 nlogn body_2000.data > nlogn_1_2000.out
# mpirun -np 7 nlogn body_2000.data > nlogn_7_2000.out
# mpirun -np 8 nlogn body_2000.data > nlogn_8_2000.out

# mpirun -np 1 nlogn body_2500.data > nlogn_1_2500.out
# mpirun -np 7 nlogn body_2500.data > nlogn_7_2500.out
# mpirun -np 8 nlogn body_2500.data > nlogn_8_2500.out