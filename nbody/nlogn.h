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

#include "cell.h"

using namespace std;

// Return current time, for performance measurement
// Adapt from previous assignment
uint64_t GetTimeStamp() {
    struct timeval tv;
    gettimeofday(&tv,NULL);
    return tv.tv_sec*(uint64_t)1000000+tv.tv_usec;
}

int N {0};
int T {0};

void load_data(string filename, vector<Body> &n_bodies)
{
    std::ifstream file {filename};
    string line {};

    //read the particals number
    std::getline(file, line);
    std::istringstream num_iss {line};

    // n bodies
    num_iss >> N;
    cout << "get " << N << " particals" << endl;

    // read the iteration values
    std::getline(file, line);
    std::istringstream iteration_iss {line};
    iteration_iss >> T;
    cout << "run " << T << " iterations" << endl;

    int count {0};

    // read each partical settings
    while (std::getline(file, line))
    {   
        std::istringstream iss {line};
        double mass, px, py, pz, vx, vy, vz;
        if (!(iss >> mass >> px >> py >> pz >> vx >> vy >> vz)) 
        {
            break; 
        }

        n_bodies.emplace_back(count, mass, px, py, pz, vx, vy, vz);
        count++;
    }

    if (count != N) {
        cout << "bodies number unmatched" << " get " << count << ", need " << N << endl;
        exit(1);
    }
}

#endif // __NLOGN_H__