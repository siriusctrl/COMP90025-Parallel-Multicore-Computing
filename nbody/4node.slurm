#!/bin/bash
#SBATCH --error=mpi-4n1c-%j.err
#SBATCH --output=mpi-4n1c-%j.out
#SBATCH --partition=snowy
#SBATCH --time=01:00:00
#SBATCH --nodes=4
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=32G #SBATCH --job-name=mpi

module load gcc/8.3.0
module load openmpi/3.1.4

mpicxx -std=c++14 -fopenmp -o nlogn nlogn.h cell.h nlogn.cpp partical.h -O3
mpirun -np 4 nlogn body_10.data > 4n1c_10_nlogn.out
mpirun -np 4 nlogn body_100.data > 4n1c_100_nlogn.out
mpirun -np 4 nlogn body_1000.data > 4n1c_1000_nlogn.out
mpirun -np 4 nlogn body_5000.data > 4n1c_5000_nlogn.out