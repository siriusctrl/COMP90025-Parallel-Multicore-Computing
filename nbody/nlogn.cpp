#include "nlogn.h"

// compute force for body i based on body j where j != i
// i: body i
void compute_force(int i, int N, Partical *particals, Force * force, Cell* cell){
    // reset force
    force->fx = 0;
    force->fy = 0;
    force->fz = 0;

    compute_force_from_octtree(cell, i, particals, force);
}

// update position and velocity of a body
void update_body(Partical * body_next, int N, Partical body_cur, Force body_force) {
    // factor = dt / m
    double factor = DT / body_cur.mass;

    body_next->px += DT * body_cur.vx;
    body_next->py += DT * body_cur.vy;
    body_next->pz += DT * body_cur.vz;

    // v = dt * a = dt * (F / m) = (dt / m) * F = factor * F
    body_next->vx += factor * body_force.fx;
    body_next->vy += factor * body_force.fy;
    body_next->vz += factor * body_force.fz;

    // wrap the velocity if out of bound
    if (body_next->px >= X_BOUND) {
        body_next->vx = -1 * abs(body_next->vx);
    } else if (body_next->px <= 0) {
        body_next->vx = abs(body_next->vx);
    }

    if (body_next->py >= Y_BOUND) {
        body_next->vy = -1 * abs(body_next->vy);
    } else if (body_next->py <= 0) {
        body_next->vy = abs(body_next->vy);
    }

    if (body_next->pz >= Z_BOUND) {
        body_next->vz = -1 * abs(body_next->vz);
    } else if (body_next->pz <= 0) {
        body_next->vz = -1 * abs(body_next->vz);
    }
}

void calculate(int N, int T, Partical *particals) {
    Partical particals_next[N];
    for (int i = 0; i < N; ++i) {
        particals_next[i] = particals[i];
    }

    Force particals_forces[N];
    for (int z = 0; z < T; ++z) {
        Cell* octtree = generate_octtree(N, particals);
        // cout << "tree generated " << octtree << endl;
        octtree->compute_cell_properties(particals);
        for (int i = 0; i < N; ++i) {
            compute_force(i, N, particals, &(particals_forces[i]), octtree);
        }

        for (int i = 0; i < N; ++i) {
            update_body(&(particals_next[i]), N, particals[i], particals_forces[i]);
        }

        for (int i = 0; i < N; i++) {
            particals[i] = particals_next[i];
        }
        delete_octtree(octtree);
    }
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    int rank;
    MPI_Comm_rank(comm, &rank);

    uint64_t start, end;

    if (argc < 2) {
        cout << "unspecified file name" << endl;
    }

    vector<Partical> particals {};
    load_data(argv[1], particals);

    start = GetTimeStamp();
    calculate(N, T, &particals.front());


    for (auto const &b: particals)
    {
        cout << b;
    }

    cout << "time = " << GetTimeStamp() - start << " ns" << endl;

    MPI::Finalize();
    return 0;
}