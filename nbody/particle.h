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


constexpr double G = 6.67e-11;
constexpr double DT = 0.001;
constexpr double MASS_BOUND = 1.0e24;

// // Particle related calculation
// class Particle {
// public:
//     int count;
//     double mass, px, py, pz, vx, vy, vz;
//     // it does not compile if we add constructor here
//     Particle() = default;
//     Particle(double init) 
//     {   
//         Particle(init, init, init, init, init, init, init, init);
//     }
//     Particle(int count, double mass, double px, double py, double pz, double vx, double vy, double vz) 
//         : count {count}, mass {mass*MASS_BOUND}, px {px * X_BOUND}, py {py * Y_BOUND}, pz {pz * Z_BOUND}, vx {vx}, vy {vy}, vz {vz}
//     {
//         // cout << "particle init successfully" << endl;
//     }

//     double compute_distance(const Particle &other) const
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

//     void set_coordinates(const Particle &other, double x_offset, double y_offset, double z_offset)
//     {
//         px = other.px + x_offset;
//         py = other.py + y_offset;
//         pz = other.pz + z_offset;
//     }

//     void set_coordinates(const Particle &other)
//     {
//         set_coordinates(other, 0.0, 0.0, 0.0);
//     }
// };

struct Particle {
    double mass, px, py, pz, vx, vy, vz;
};

Particle default_particle()
{
    return Particle {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
}

Particle init_particle(double mass, double px, double py, double pz, double vx, double vy, double vz) {
    return Particle {mass * MASS_BOUND, px * X_BOUND, py * Y_BOUND, pz * Z_BOUND, vx, vy, vz};
}

double compute_distance(const Particle &a, const Particle &b) {
        double dx, dy, dz;
        dx = a.px - b.px;
        dy = a.py - b.py;
        dz = a.pz - b.pz;
        return sqrt(pow(dx, 2) + pow(dy, 2) + pow(dz, 2));
}

void set_coordinates(Particle &target, const Particle &other, double x_offset, double y_offset, double z_offset)
{
    target.px = other.px + x_offset;
    target.py = other.py + y_offset;
    target.pz = other.pz + z_offset;
}

void set_coordinates(Particle &target, const Particle &other)
{
    set_coordinates(target, other, 0.0, 0.0, 0.0);
}

// Overloaded operator for '<<' for struct output
std::ostream& operator<<(std::ostream &strm, const Particle &particle) {
    strm << "mass: " << particle.mass << " px: " << particle.px << " py: " << particle.py << " pz: " << particle.pz;
    strm << " vx: " << particle.vx << " vy: " << particle.vy << " vz: " << particle.vz << endl;
    return strm;
}

class Force {
public:
    double fx, fy, fz;
};
#endif // __PARTICAL_H__