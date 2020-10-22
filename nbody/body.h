#ifndef __BODY_H__
#define __BODY_H__

#include <cmath>
#include <string>
#include <sstream>

using namespace std;

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
#endif // __BODY_H__