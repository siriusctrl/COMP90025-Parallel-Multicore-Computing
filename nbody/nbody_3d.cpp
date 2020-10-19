#include "nbody_3d.h"

int ITERATION = 0;

int main(int argc, char const *argv[])
{
    if (argc < 2)
    {
        cout << "please specify data file" << endl;
        return 0;
    }

    vector<partical> particals;
    read_data(argv[1], particals);

    omp_set_num_threads(12);

    for (int t=0; t<ITERATION; t++)
    {
        vector<partical> next_particals {particals};

        #pragma omp parallel for
        for (int i=0; i < particals.size(); ++i)
        {
            update_coordinate(particals, i, next_particals);
            // cout << "\nThis:\n" << endl;
            // for (auto const &p:particals) 
            // {
            //     cout << p;
            // }

            // cout << "\nNest:\n" << endl;
            // for (auto const &p:next_particals) 
            // {
            //     cout << p;
            // }
        }

        particals = next_particals;

        // cout << "----------------------------------" << endl;
    }

    cout << "\nEnd:\n" << endl;

    for (auto const &p:particals) 
    {
        cout << p;
    }

    return 0;
}

void update_coordinate(const vector<partical> &particals, 
                        int index, 
                        vector<partical> &next_particals)
{
    double epsilon {0.001};
    double dx {0.0}, dy {0.0}, dz {0.0};
    double r {0.0}, r_square {0.0}, f_share {0.0}, fx {0.0}, fy{0.0}, fz {0.0};

    for (int i=0; i < particals.size(); ++i)
    {
        if (i != index)
        {
            dx = particals.at(index).px - particals.at(i).px;
            dy = particals.at(index).py - particals.at(i).py;
            dz = particals.at(index).pz - particals.at(i).pz;

            r_square = std::sqrt(std::pow(dx,2)) + 
                        std::sqrt(std::pow(dy, 2)) + 
                        std::sqrt(std::pow(dz, 2)) + epsilon;

            r = std::sqrt(r_square) + epsilon;

            f_share = G * particals.at(i).mass * particals.at(index).mass / r_square;

            fx = f_share * (dx) / r;
            fy = f_share * (dy) / r;
            fz = f_share * (dz) / r;

            next_particals.at(index).px += DT * particals.at(index).px;
            next_particals.at(index).py += DT * particals.at(index).py;
            next_particals.at(index).pz += DT * particals.at(index).pz;

            next_particals.at(index).vx += fx * DT / particals.at(index).mass;
            next_particals.at(index).vy += fy * DT / particals.at(index).mass;
            next_particals.at(index).vz += fz * DT / particals.at(index).mass;
        }
    }
}

void read_data(string filename, vector<partical> &particals) 
{

    std::ifstream file {filename};
    string line {};

    int num {0};

    //read the particals number
    std::getline(file, line);
    std::istringstream num_iss {line};
    num_iss >> num;
    cout << "get " << num << " particals" << endl;
    
    // read the iteration
    std::getline(file, line);
    std::istringstream iteration_iss {line};
    iteration_iss >> ITERATION;
    cout << "run " << ITERATION << " iterations" << endl;

    double mass, px, py, pz, vx, vy, vz;

    // read each partical settings
    while (std::getline(file, line))
    {   
        std::istringstream iss {line};
        double mass, px, py, pz, vx, vy, vz;
        if (!(iss >> mass >> px >> py >> pz >> vx >> vy >> vz)) 
        {
            break; 
        }

        particals.emplace_back(mass, px, py, pz, vx, vy, vz);
    }

    if (particals.size() != num) {
        cout << "partical number unmatched" << endl;
        exit(1);
    }
}

// g++ -std=c++14 -fopenmp nbody_3d.cpp -o nbody_3d -O3