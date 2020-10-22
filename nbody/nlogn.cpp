#include "nlogn.h"

// compute force for body i based on body j where j != i
// i: body i
void compute_force(int i, int N, double G, Body *n_bodies, Force * force, Cell* cell){
    // reset force
    force->fx = 0;
    force->fy = 0;
    force->fz = 0;

    compute_force_from_octtree(cell, i, n_bodies, G, force);
}

// update position and velocity of a body
void update_body(Body * body_next, int N, double G, double TIME_DELTA, Body body_cur, Force body_force) {
    // factor = dt / m
    double factor = TIME_DELTA / body_cur.mass;

    body_next->px += TIME_DELTA * body_cur.vx;
    body_next->py += TIME_DELTA * body_cur.vy;
    body_next->pz += TIME_DELTA * body_cur.vz;

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

void calculate(int N, int T, double G, double TIME_DELTA, Body *n_bodies) {
    Body n_bodies_next[N];
    for (int i = 0; i < N; ++i) {
        n_bodies_next[i] = n_bodies[i];
    }

    Force n_bodies_forces[N];
    for (int z = 0; z < T; ++z) {
        Cell* octree = generate_octtree(N, n_bodies);
        // cout << "tree generated " << octree << endl;
        compute_cell_properties(octree, n_bodies);
        for (int i = 0; i < N; ++i) {
            compute_force(i, N, G, n_bodies, &(n_bodies_forces[i]), octree);
        }
        // cout << "force computed" << endl;
        for (int i = 0; i < N; ++i) {
            update_body(&(n_bodies_next[i]), N, G, TIME_DELTA, n_bodies[i], n_bodies_forces[i]);
        }

        for (int i = 0; i < N; i++) {
            n_bodies[i] = n_bodies_next[i];
        }
        delete_octtree(octree);
    }
}

// void calculate(int N, int T, double G, double TIME_DELTA, vector<Body> n_bodies) {
//     Body n_bodies_next[N];
//     // vector<Body> n_bodies_next {n_bodies}
//     for (int i = 0; i < N; ++i) {
//         n_bodies_next[i] = (&n_bodies.front())[i];
//     }

//     Force n_bodies_forces[N];
//     for (int z = 0; z < T; ++z) {
//         Cell* octree = generate_octtree(N, &n_bodies.front());
//         // cout << "tree generated " << octree << endl;
//         compute_cell_properties(octree, &n_bodies.front());

//         for (int i = 0; i < N; ++i) {
//             compute_force(i, N, G, &n_bodies.front(), &(n_bodies_forces[i]), octree);
//         }
//         // cout << "force computed" << endl;
//         for (int i = 0; i < N; ++i) {
//             update_body(&(n_bodies_next[i]), N, G, TIME_DELTA, n_bodies[i], n_bodies_forces[i]);
//         }

//         for (int i = 0; i < N; i++) {
//             n_bodies[i] = n_bodies_next[i];
//         }

//         delete_octtree(octree);
//     }
// }

int main(int argc, char **argv) {
    uint64_t start, end;

    if (argc < 2) {
        cout << "unspecified file name" << endl;
    }

    vector<Body> n_bodies {};
    load_data(argv[1], n_bodies);

    start = GetTimeStamp();
    calculate(N, T, G, DT, &n_bodies.front());
    // calculate(N, T, G, DT, n_bodies);


    for (auto const &b: n_bodies) 
    {
        cout << b;
    }

    cout << "time = " << GetTimeStamp() - start << " ns" << endl;

    return 0;
}