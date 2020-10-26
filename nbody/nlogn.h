#ifndef __NLOGN_H__
#define __NLOGN_H__

#include <stdio.h>
#include <sys/time.h>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <vector>
#include <sstream>
#include <fstream>
#include <mpi.h>
#include <omp.h>

#include "cell.h"

using std::string;
using std::vector;

// Return current time, for performance measurement
// Adapt from previous assignment
uint64_t GetTimeStamp() {
    struct timeval tv;
    gettimeofday(&tv,NULL);
    return tv.tv_sec*(uint64_t)1000000+tv.tv_usec;
}

int N {0};
int T {0};
const MPI_Comm comm = MPI_COMM_WORLD;
constexpr int root {0};

void load_data(string filename, vector<Particle> &particles)
{
    std::ifstream file {filename};
    string line {};

    //read the particles number
    std::getline(file, line);
    std::istringstream num_iss {line};

    // n bodies
    num_iss >> N;
    std::cout << "get " << N << " particles" << endl;

    // read the iteration values
    std::getline(file, line);
    std::istringstream iteration_iss {line};
    iteration_iss >> T;
    std::cout << "run " << T << " iterations" << endl;

    int count {0};

    // read each particle settings
    while (std::getline(file, line))
    {   
        std::istringstream iss {line};
        double mass, px, py, pz, vx, vy, vz;
        if (!(iss >> mass >> px >> py >> pz >> vx >> vy >> vz)) 
        {
            break; 
        }

        particles.emplace_back(init_particle(mass, px, py, pz, vx, vy, vz));
        count++;
    }

    if (count != N) {
        std::cout << "bodies number unmatched" << " get " << count << ", need " << N << endl;
        exit(1);
    }
}

MPI_Datatype create_MPI_Particle() {
    MPI_Datatype MPI_Particle;
    MPI_Type_contiguous(7, MPI_DOUBLE, &MPI_Particle);
    MPI_Type_commit(&MPI_Particle);
    return MPI_Particle;
}

MPI_Datatype create_MPI_Force() {
    MPI_Datatype MPI_Force;
    MPI_Type_contiguous(3, MPI_DOUBLE, &MPI_Force);
    MPI_Type_commit(&MPI_Force);
    return MPI_Force;
}

#endif // __NLOGN_H__