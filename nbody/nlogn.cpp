#include "nlogn.h"

/*
 * Compute force to every other particals, excluding itself.
 */
void compute_force(int i, int N, Partical *particals, Force * force, Cell* cell){
    // reset force
    force->fx = 0;
    force->fy = 0;
    force->fz = 0;

    BH_Octtree::octtree_force(cell, i, particals, force);
}

// update position and velocity of a body
void update_particals(Partical *next_partical, int N, const Partical &cur_partical, Force force) {
    // factor = dt / m
    double factor = DT / cur_partical.mass;

    next_partical->px += DT * cur_partical.vx;
    next_partical->py += DT * cur_partical.vy;
    next_partical->pz += DT * cur_partical.vz;

    // v = dt * a = dt * (F / m) = (dt / m) * F = factor * F
    next_partical->vx += factor * force.fx;
    next_partical->vy += factor * force.fy;
    next_partical->vz += factor * force.fz;

    // wrap the velocity if out of bound to simplify the problem
    if (next_partical->px >= X_BOUND) {
    next_partical->vx = -1 * abs(next_partical->vx);
    } else if (next_partical->px <= 0) {
    next_partical->vx = abs(next_partical->vx);
    }

    if (next_partical->py >= Y_BOUND) {
    next_partical->vy = -1 * abs(next_partical->vy);
    } else if (next_partical->py <= 0) {
    next_partical->vy = abs(next_partical->vy);
    }

    if (next_partical->pz >= Z_BOUND) {
    next_partical->vz = -1 * abs(next_partical->vz);
    } else if (next_partical->pz <= 0) {
    next_partical->vz = -1 * abs(next_partical->vz);
    }
}

void calculate(Partical *particals) {
    MPI_Datatype MPI_Partical = create_MPI_Partical();
    MPI_Datatype MPI_Force = create_MPI_Force();
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    int workload = (int) ceil((1.0 * N) / size);
    int start    = rank * workload;
    int end      = std::min(N, (rank+1) * workload);

    Force total_forces[workload * size];
    Partical total_particals[workload * size];

    if (rank == root) {
        for (int i=0; i < N; ++i) {
            total_particals[i] = particals[i];
        }
    }

    MPI_Bcast(total_particals, N, MPI_Partical, root, comm);

    // omp_set_num_threads(1);

    //repeat the following for T iterations
    for (int adc = 0; adc < T; ++adc) 
    {
        if(rank == root) {
            std::cout << adc << endl;
        }
        
        Force current_forces[workload];
        Partical current_particals[workload];
        Cell* octtree = BH_Octtree::create_tree(N, total_particals);
        octtree->compute_cell_properties(total_particals);

        // #pragma omp parallel for
        for (int i = start; i < end; ++i) 
        {
            current_particals[i - start] = total_particals[i];
        }

        // #pragma omp parallel for
        for (int i = start; i < end; ++i)
        {
            compute_force(i, N, total_particals, current_forces+i-start, octtree);
        }

        MPI_Allgather(current_forces, workload, MPI_Force,
                        total_forces, workload, MPI_Force,
                        comm);
        
        // #pragma omp parallel for
        for (int i = start; i < end; ++i)
        {
            update_particals(current_particals+i-start, N, total_particals[i], total_forces[i]);
        }

        MPI_Allgather(current_particals, workload, MPI_Partical,
                        total_particals, workload, MPI_Partical,
                        comm);

        BH_Octtree::delete_octtree(octtree);
    }

    MPI_Type_free(&MPI_Partical);
    MPI_Type_free(&MPI_Force);
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    
    int rank;
    MPI_Comm_rank(comm, &rank);

    uint64_t start, end;
    vector<Partical> particals {};

    if (rank == 0 ) 
    {
        if (argc < 2) {
            std::cout << "unspecified file name" << endl;
        }


        load_data(argv[1], particals);

        start = GetTimeStamp();
    }

    MPI_Bcast(&N, 1, MPI_INT, root, comm);
    MPI_Bcast(&T, 1, MPI_INT, root, comm);
    calculate(&particals.front());

    if (rank == 0) {
        for (auto const &b: particals)
        {
            std::cout << b;
        }

        std::cout << "time = " << (GetTimeStamp() - start) / 1e6 << " s" << endl;
    }

    MPI_Finalize();
    return 0;
}
//mpicxx -std=c++14 -fopenmp -o nlogn nlogn.h cell.h nlogn.cpp partical.h -O3  