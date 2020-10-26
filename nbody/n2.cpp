#include "particle.h"
#include "n2.h"

int ITERATION = 0;

int main(int argc, char const *argv[])
{
    if (argc < 2)
    {
        cout << "please specify data file" << endl;
        return 0;
    }

    vector<particle> particles;
    read_data(argv[1], particles);

    omp_set_num_threads(12);

    for (int t=0; t<ITERATION; t++)
    {
        vector<particle> next_particles {particles};

        #pragma omp parallel for
        for (int i=0; i < particles.size(); ++i)
        {
            update_coordinate(particles, i, next_particles);
            // cout << "\nThis:\n" << endl;
            // for (auto const &p:particles) 
            // {
            //     cout << p;
            // }

            // cout << "\nNest:\n" << endl;
            // for (auto const &p:next_particles) 
            // {
            //     cout << p;
            // }
        }

        particles = next_particles;

        // cout << "----------------------------------" << endl;
    }

    cout << "\nEnd:\n" << endl;

    for (auto const &p:particles) 
    {
        cout << p;
    }

    return 0;
}

void update_coordinate(const vector<particle> &particles, 
                        int index, 
                        vector<particle> &next_particles)
{
    double epsilon {0.001};
    double dx {0.0}, dy {0.0}, dz {0.0};
    double r {0.0}, r_square {0.0}, f_share {0.0}, fx {0.0}, fy{0.0}, fz {0.0};

    for (int i=0; i < particles.size(); ++i)
    {
        if (i != index)
        {
            dx = particles.at(index).px - particles.at(i).px;
            dy = particles.at(index).py - particles.at(i).py;
            dz = particles.at(index).pz - particles.at(i).pz;

            r_square = std::sqrt(std::pow(dx,2)) + 
                        std::sqrt(std::pow(dy, 2)) + 
                        std::sqrt(std::pow(dz, 2)) + epsilon;

            r = std::sqrt(r_square) + epsilon;

            f_share = G * particles.at(i).mass * particles.at(index).mass / r_square;

            fx = f_share * (dx) / r;
            fy = f_share * (dy) / r;
            fz = f_share * (dz) / r;

            next_particles.at(index).px += DT * particles.at(index).px;
            next_particles.at(index).py += DT * particles.at(index).py;
            next_particles.at(index).pz += DT * particles.at(index).pz;

            next_particles.at(index).vx += fx * DT / particles.at(index).mass;
            next_particles.at(index).vy += fy * DT / particles.at(index).mass;
            next_particles.at(index).vz += fz * DT / particles.at(index).mass;

            particle* next_particle = &next_particles.at(index);

            // wrap the velocity if out of bound to simplify the problem
            if (next_particle->px >= X_BOUND) {
            next_particle->vx = -1 * abs(next_particle->vx);
            } else if (next_particle->px <= 0) {
            next_particle->vx = abs(next_particle->vx);
            }

            if (next_particle->py >= Y_BOUND) {
            next_particle->vy = -1 * abs(next_particle->vy);
            } else if (next_particle->py <= 0) {
            next_particle->vy = abs(next_particle->vy);
            }

            if (next_particle->pz >= Z_BOUND) {
            next_particle->vz = -1 * abs(next_particle->vz);
            } else if (next_particle->pz <= 0) {
            next_particle->vz = -1 * abs(next_particle->vz);
            }
        }
    }
}

void read_data(string filename, vector<Particle> &particles) 
{

    std::ifstream file {filename};
    string line {};

    int num {0};

    //read the particles number
    std::getline(file, line);
    std::istringstream num_iss {line};
    num_iss >> num;
    cout << "get " << num << " particles" << endl;
    
    // read the iteration
    std::getline(file, line);
    std::istringstream iteration_iss {line};
    iteration_iss >> ITERATION;
    cout << "run " << ITERATION << " iterations" << endl;

    double mass, px, py, pz, vx, vy, vz;

    // read each particle settings
    while (std::getline(file, line))
    {   
        std::istringstream iss {line};
        double mass, px, py, pz, vx, vy, vz;
        if (!(iss >> mass >> px >> py >> pz >> vx >> vy >> vz)) 
        {
            break; 
        }

        particles.emplace_back(mass, px, py, pz, vx, vy, vz);
    }

    if (particles.size() != num) {
        cout << "particle number unmatched" << endl;
        exit(1);
    }
}

// g++ -std=c++14 -fopenmp nbody_3d.cpp -o nbody_3d -O3