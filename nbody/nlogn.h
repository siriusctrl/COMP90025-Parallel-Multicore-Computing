#ifndef __NLOGN_H__
#define __NLOGN_H__

#include <stdio.h>
#include <sys/time.h>
#include <math.h>
#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <vector>
#include <sstream>
#include <fstream>

using namespace std;

// Return current time, for performance measurement
// Adapt from previous assignment
uint64_t GetTimeStamp() {
    struct timeval tv;
    gettimeofday(&tv,NULL);
    return tv.tv_sec*(uint64_t)1000000+tv.tv_usec;
}

constexpr double EPSILON = 0.000001;
constexpr double X_BOUND = 1.0e6;      // Width of space
constexpr double Y_BOUND = 1.0e6;      // Height of space
constexpr double Z_BOUND = 1.0e6;      // Depth of space
constexpr double THETA   = 0;        // Opening angle, for approximation in Barned hut algorithm

constexpr double G = 6.67e-11;
constexpr double DT = 0.001;
constexpr double MASS_BOUND = 1.0e24;

int N {0};
int T {0};

// Body related calculation
class Body {
public:
    int count;
    double mass, px, py, pz, vx, vy, vz;
    // it does not compile if we add constructor here
    Body() = default;
    Body(double init) 
    {   
        Body(init, init, init, init, init, init, init, init);
    }
    Body(int count, double mass, double px, double py, double pz, double vx, double vy, double vz) 
        : count {count}, mass {mass*MASS_BOUND}, px {px * X_BOUND}, py {py * Y_BOUND}, pz {pz * Z_BOUND}, vx {vx}, vy {vy}, vz {vz}
    {
        // cout << "body init successfully" << endl;
    }

    double compute_distance(const Body &other)
    {
        double dx, dy, dz;
        dx = other.px - px;
        dy = other.py - py;
        dz = other.pz - pz;
        return sqrt(pow(dx, 2) + pow(dy, 2) + pow(dz, 2));
    }

    void set_coordinates(double px, double py, double pz)
    {
        this->px = px;
        this->py = py;
        this->pz = pz;
    }

    void set_coordinates(const Body &other, double x_offset, double y_offset, double z_offset)
    {
        px = other.px + x_offset;
        py = other.py + y_offset;
        pz = other.pz + z_offset;
    }

    void set_coordinates(const Body &other)
    {
        set_coordinates(other, 0.0, 0.0, 0.0);
    }
};

// Overloaded operator for '<<' for struct output
std::ostream& operator<<(std::ostream &strm, const Body &body) {
    strm << body.count << ": ";
    strm << "mass: " << body.mass << " px: " << body.px << " py: " << body.py << " pz: " << body.pz;
    strm << " vx: " << body.vx << " vy: " << body.vy << " vz: " << body.vz << endl;
    return strm;
}

class Force {
public:
    double fx, fy, fz;
};

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