#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <omp.h>
#include <cmath>

using std::vector;
using std::string;
using std::cout;
using std::endl;

constexpr double X_BOUND = 1.0e6;      // Width of space
constexpr double Y_BOUND = 1.0e6;      // Height of space
constexpr double Z_BOUND = 1.0e6;      // Depth of space
constexpr double G = 6.67e-11;
constexpr double DT = 0.001;
constexpr double MASS_BOUND = 1.0e24;

class particle
{
public:
    double mass, px, py, pz, vx, vy, vz;
    particle(double mass, double px, double py, double pz, double vx, double vy, double vz) 
            : mass {mass*MASS_BOUND}, px {px * X_BOUND}, py {py * Y_BOUND}, pz {pz * Z_BOUND}, vx {vx}, vy {vy}, vz {vz}
    {
        // cout << "particle init successfully" << endl;
    }
};

// easy print the internal setting for class particle
std::ostream& operator<<(std::ostream &strm, const particle &p) {
    strm << "mass: " << p.mass << " px: " << p.px << " py: " << p.py << " pz: " << p.pz;
    strm << " vx: " << p.vx << " vy: " << p.vy << " vz: " << p.vz << endl;
    return strm;
}

void update_coordinate(const vector<particle> &particles, int index, vector<particle> &next_particles);
// read particle config from given config file
void read_data(string filename, vector<particle> &particles);