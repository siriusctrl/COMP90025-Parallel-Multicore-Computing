#include "n2.h"

int ITERATION = 0;

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    
    int rank;
    MPI_Comm_rank(comm, &rank);

    uint64_t start, end;
    vector<Particle> particles {};

    if (rank == 0 ) 
    {
        if (argc < 2) {
            std::cout << "unspecified file name" << endl;
        }


        load_data(argv[1], particles);

        start = GetTimeStamp();
    }

    MPI_Bcast(&N, 1, MPI_INT, root, comm);
    MPI_Bcast(&T, 1, MPI_INT, root, comm);
    simulate(&particles.front());

    if (rank == 0) {
        for (auto const &b: particles)
        {
            std::cout << b;
        }

        std::cout << "time = " << (GetTimeStamp() - start) / 1e6 << " s" << endl;
    }

    MPI_Finalize();
    return 0;
}

void simulate(Particle *particles) {
    MPI_Datatype MPI_Particle = create_MPI_Particle();
    MPI_Datatype MPI_Force = create_MPI_Force();
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    int workload = (int) ceil((1.0 * N) / size);
    int start    = rank * workload;
    int end      = std::min(N, (rank+1) * workload);

    Force total_forces[workload * size];
    Particle total_particles[workload * size];

    if (rank == root) {
        for (int i=0; i < N; ++i) {
            total_particles[i] = particles[i];
        }
    }

    MPI_Bcast(total_particles, N, MPI_Particle, root, comm);

    // omp_set_num_threads(1);

    //repeat the following for T iterations
    for (int iter = 0; iter < T; ++iter) 
    {   

        Force current_forces[workload];
        Particle current_particles[workload];

        #pragma omp parallel for num_threads(2)
        for (int i = start; i < end; ++i) 
        {
            current_particles[i - start] = total_particles[i];
        }

        #pragma omp parallel for num_threads(2)
        for (int i = start; i < end; ++i)
        {
            compute_force(i, N, total_particles, current_forces+i-start);
        }

        MPI_Allgather(current_forces, workload, MPI_Force,
                        total_forces, workload, MPI_Force,
                        comm);
        
        // #pragma omp parallel for
        for (int i = start; i < end; ++i)
        {
            update_particles(current_particles+i-start, total_particles[i], total_forces[i]);
        }

        MPI_Allgather(current_particles, workload, MPI_Particle,
                        total_particles, workload, MPI_Particle,
                        comm);
    }

    MPI_Type_free(&MPI_Particle);
    MPI_Type_free(&MPI_Force);
}

void compute_force(int i, int N, Particle *particles, Force * force) 
{
    force->fx = 0;
    force->fy = 0;
    force->fz = 0;

    for (int p=0; p < N ; ++p)
    {
        if (p != i) 
        {
            double px_diff, py_diff, pz_diff, factor, d;
            // distance in x direction
            px_diff = particles[p].px - particles[i].px;
            // distance in y direction
            py_diff = particles[p].py - particles[i].py;
            // distance in z direction
            pz_diff = particles[p].pz - particles[i].pz;

            d = compute_distance(particles[i], particles[p]);

            // even though we are aussing that no collide, but we still add a small number to
            // prevent 0 division
            factor = G * particles[p].mass * particles[i].mass / (pow(d, 3) + EPSILON);
            force->fx += px_diff * factor; // force in x direction
            force->fy += py_diff * factor; // force in y direction
            force->fz += pz_diff * factor; // force in z direction 
        }
    }
}

void update_particles(Particle *next_particle, const Particle &cur_particle, Force force) 
{
    // factor = dt / m
    double factor = DT / cur_particle.mass;

    next_particle->px += DT * cur_particle.vx;
    next_particle->py += DT * cur_particle.vy;
    next_particle->pz += DT * cur_particle.vz;

    // v = dt * a = dt * (F / m) = (dt / m) * F = factor * F
    next_particle->vx += factor * force.fx;
    next_particle->vy += factor * force.fy;
    next_particle->vz += factor * force.fz;

    // wrap the velocity if out of "bound"
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

// mpicxx -std=c++14 -fopenmp n2.cpp n2.h particle.h -o n2 -O3