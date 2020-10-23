#ifndef __PARTICAL_H__
#define __PARTICAL_H__

#include <cmath>
#include <string>
#include <sstream>

using std::endl;

constexpr double EPSILON = 0.000001;
constexpr double X_BOUND = 1.0e6;      // Width of space
constexpr double Y_BOUND = 1.0e6;      // Height of space
constexpr double Z_BOUND = 1.0e6;      // Depth of space
constexpr double THETA   = 0;        // Opening angle, for approximation in Barned hut algorithm

constexpr double G = 6.67e-11;
constexpr double DT = 0.001;
constexpr double MASS_BOUND = 1.0e24;

// // Partical related calculation
// class Partical {
// public:
//     int count;
//     double mass, px, py, pz, vx, vy, vz;
//     // it does not compile if we add constructor here
//     Partical() = default;
//     Partical(double init) 
//     {   
//         Partical(init, init, init, init, init, init, init, init);
//     }
//     Partical(int count, double mass, double px, double py, double pz, double vx, double vy, double vz) 
//         : count {count}, mass {mass*MASS_BOUND}, px {px * X_BOUND}, py {py * Y_BOUND}, pz {pz * Z_BOUND}, vx {vx}, vy {vy}, vz {vz}
//     {
//         // cout << "partical init successfully" << endl;
//     }

//     double compute_distance(const Partical &other) const
//     {
//         double dx, dy, dz;
//         dx = other.px - px;
//         dy = other.py - py;
//         dz = other.pz - pz;
//         return sqrt(pow(dx, 2) + pow(dy, 2) + pow(dz, 2));
//     }

//     void set_coordinates(double px, double py, double pz)
//     {
//         this->px = px;
//         this->py = py;
//         this->pz = pz;
//     }

//     void set_coordinates(const Partical &other, double x_offset, double y_offset, double z_offset)
//     {
//         px = other.px + x_offset;
//         py = other.py + y_offset;
//         pz = other.pz + z_offset;
//     }

//     void set_coordinates(const Partical &other)
//     {
//         set_coordinates(other, 0.0, 0.0, 0.0);
//     }
// };

struct Partical {
    double mass, px, py, pz, vx, vy, vz;
};

Partical default_partical()
{
    return Partical {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
}

Partical init_partical(double mass, double px, double py, double pz, double vx, double vy, double vz) {
    return Partical {mass * MASS_BOUND, px * X_BOUND, py * Y_BOUND, pz * Z_BOUND, vx, vy, vz};
}

double compute_distance(const Partical &a, const Partical &b) {
        double dx, dy, dz;
        dx = a.px - b.px;
        dy = a.py - b.py;
        dz = a.pz - b.pz;
        return sqrt(pow(dx, 2) + pow(dy, 2) + pow(dz, 2));
}

void set_coordinates(Partical &target, const Partical &other, double x_offset, double y_offset, double z_offset)
{
    target.px = other.px + x_offset;
    target.py = other.py + y_offset;
    target.pz = other.pz + z_offset;
}

void set_coordinates(Partical &target, const Partical &other)
{
    set_coordinates(target, other, 0.0, 0.0, 0.0);
}

// Overloaded operator for '<<' for struct output
std::ostream& operator<<(std::ostream &strm, const Partical &partical) {
    strm << "mass: " << partical.mass << " px: " << partical.px << " py: " << partical.py << " pz: " << partical.pz;
    strm << " vx: " << partical.vx << " vy: " << partical.vy << " vz: " << partical.vz << endl;
    return strm;
}

class Force {
public:
    double fx, fy, fz;
};
#endif // __PARTICAL_H__