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

// Body related calculation
class Body {
public:
    double mass, px, py, pz, vx, vy, vz;
    // it does not compile if we add constructor here
    Body() = default;
    Body(double mass, double px, double py, double pz, double vx, double vy, double vz) 
            : mass {mass*MASS_BOUND}, px {px * X_BOUND}, py {py * Y_BOUND}, pz {pz * Z_BOUND}, vx {vx}, vy {vy}, vz {vz}
    {
        // cout << "body init successfully" << endl;
    }
};

// Overloaded operator for '<<' for struct output
std::ostream& operator<<(std::ostream &strm, const Body &body) {
    strm << "mass: " << body.mass << " px: " << body.px << " py: " << body.py << " pz: " << body.pz;
    strm << " vx: " << body.vx << " vy: " << body.vy << " vz: " << body.vz;
    return strm;
}

struct Force {
    double fx, fy, fz;
};

#endif // __NLOGN_H__