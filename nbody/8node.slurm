#!/bin/bash
#SBATCH --error=mpi-8n1c-%j.err
#SBATCH --output=mpi-8n1c-%j.extra
#SBATCH --partition=snowy
#SBATCH --time=01:00:00
#SBATCH --nodes=8
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=32G #SBATCH --job-name=mpi

module load gcc/8.3.0
module load openmpi/3.1.4

mpicxx -std=c++14 -fopenmp -o nlogn nlogn.h cell.h nlogn.cpp partical.h -O3
# mpirun -np 8 nlogn body_10.data > 8n1c_10_nlogn.out
# mpirun -np 8 nlogn body_100.data > 8n1c_100_nlogn.out
# mpirun -np 8 nlogn body_1000.data > 8n1c_1000_nlogn.out
# mpirun -np 8 nlogn body_1500.data > 8n1c_1500_nlogn.out
# mpirun -np 8 nlogn body_2500.data > 8n1c_2500_nlogn.out
mpirun -np 8 nlogn body_5000.data > 8n1c_5000_nlogn.out
