#include <stdio.h>
#include <sys/time.h>
#include <math.h>
#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <vector>

using namespace std;

// Return current time, for performance measurement
// Adapt from previous assignment
uint64_t GetTimeStamp() {
    struct timeval tv;
    gettimeofday(&tv,NULL);
    return tv.tv_sec*(uint64_t)1000000+tv.tv_usec;
}

constexpr double EPSILON = 0.000001;
constexpr double X_BOUND = 1.0e6;      // Width of space
constexpr double Y_BOUND = 1.0e6;      // Height of space
constexpr double Z_BOUND = 1.0e6;      // Depth of space
constexpr double G = 6.67e-11;
constexpr double DT = 0.001;
constexpr double MASS_BOUND = 1.0e24;
constexpr double THETA   = 0;        // Opening angle, for approximation in Barned hut algorithm

// Body related calculation
struct Body {
    double mass, px, py, pz, vx, vy, vz;
};
// Overloaded operator for '<<' for struct output
ostream& operator << (ostream& os, const Body& body) {
    os << "body[" << body.px << "][" << body.py << "][" << body.pz << "] velocity: (" << body.vx << ", " << body.vy << ", " << body.vz << ") mass: " << body.mass;
    return os;
}

struct Force {
    double fx, fy, fz;
};

// **************************************************************************
// octtree data structure and helper methods starts
// **************************************************************************
/* Cubic cell representing tree node in Barnes-Hut algorithm */
typedef struct Cell  {
   int index;                   // Index into arrays to identify particle's 
                                // position and mass
   int n_children;              // Indicate whether cell is leaf or has 8 children
   Body center;                 // center of approximate body
                                //  - Mass of particle of total mass of subtree
                                //  - position of cell(cube) in space
                                //  - position of center of mass of cell
   double width, height, depth; // Width, Height, and Depth of cell
   struct Cell* children[8];    // Pointers to child nodes
} Cell;

/* Creates a cell to be used in the octtree */
Cell* create_cell(double width, double height, double depth) {
   Cell* cell = (Cell*) malloc(sizeof(Cell));
   cell->index = -1;
   cell->n_children = 0;
   
   (cell->center).mass = 0;
   (cell->center).px   = 0;
   (cell->center).py   = 0;
   (cell->center).pz   = 0;
   (cell->center).vx   = 0;
   (cell->center).vy   = 0;
   (cell->center).vz   = 0;

   cell->width = width;
   cell->height = height;
   cell->depth = depth;   
   
   return cell;
}

/* sets the location of the children relative to the current cell */
void set_location_of_children(Cell* cell, double width, double heigth, double depth){
   // Set location of new cells
   ((cell->children[0])->center).px = (cell->center).px;
   ((cell->children[0])->center).py = (cell->center).py;
   ((cell->children[0])->center).pz = (cell->center).pz;

   ((cell->children[1])->center).px = (cell->center).px + width;
   ((cell->children[1])->center).py = (cell->center).py;
   ((cell->children[1])->center).pz = (cell->center).pz;

   ((cell->children[2])->center).px = (cell->center).px + width;
   ((cell->children[2])->center).py = (cell->center).py;
   ((cell->children[2])->center).pz = (cell->center).pz + depth;

   ((cell->children[3])->center).px = (cell->center).px;
   ((cell->children[3])->center).py = (cell->center).py;
   ((cell->children[3])->center).pz = (cell->center).pz + depth;

   ((cell->children[4])->center).px = (cell->center).px;
   ((cell->children[4])->center).py = (cell->center).py + heigth;
   ((cell->children[4])->center).pz = (cell->center).pz;

   ((cell->children[5])->center).px = (cell->center).px + width;
   ((cell->children[5])->center).py = (cell->center).py + heigth;
   ((cell->children[5])->center).pz = (cell->center).pz;

   ((cell->children[6])->center).px = (cell->center).px + width;   // Coordinates of this cell marks
   ((cell->children[6])->center).py = (cell->center).py + heigth;  // the mid-point of the parent cell
   ((cell->children[6])->center).pz = (cell->center).pz + depth;   //
   
   ((cell->children[7])->center).px = (cell->center).px;
   ((cell->children[7])->center).py = (cell->center).py + heigth;
   ((cell->children[7])->center).pz = (cell->center).pz + depth;
}

/*
 * Generates new children for the current cell, forming a subtree. 
 * The current cell will no longer be a leaf
 */
void generate_children(Cell* cell) {
   // Calculate subcell dimensions
   double width  = cell->width / 2.0;
   double height = cell->height / 2.0;
   double depth  = cell->depth / 2.0;

   // Cell no longer a leaf
   cell->n_children = 8;   
   
   // Create and initialize new children   
   for (int i = 0; i < cell->n_children; ++i) {
      cell->children[i] = create_cell(width, height, depth);
   }
   
   set_location_of_children(cell, width, height, depth);   
}

/* Locates the child to which the particle must be added */
int locate_child(Cell* cell, Body body) {
    // Determine which child to add the body to
    if (body.px > (cell->children[6])->center.px) {
        if (body.py > (cell->children[6])->center.py) {
            if (body.pz > (cell->children[6])->center.pz) {
                return 6;
            } else {
                return 5;
            }
        } else {
            if (body.pz > (cell->children[6])->center.pz) {
                return 2;
            } else {
                return 1;
            }
        }
    } else {
        if (body.py > (cell->children[6])->center.py) {
            if (body.pz > (cell->children[6])->center.pz) {
                return 7;
            } else {
                return 4;
            }
        } else {
            if (body.pz > (cell->children[6])->center.pz) {
                return 3;
            } else {
                return 0;
            }
        }
    }
}

/* Added a particle to the cell. If a particle already
 * exists, the cube/cell is sub-divided adding the existing
 * and new particle to the sub cells
 */
void add_to_cell(Cell* cell, Body* n_bodies, int i) {
    if (cell->index == -1) {         
        cell->index = i;
        return;         
    }
         
   generate_children(cell);

   // The current cell's body must now be re-added to one of its children
   int sc1 = locate_child(cell, n_bodies[cell->index]);
   cell->children[sc1]->index = cell->index;   

   // Locate child for new body
   int sc2 = locate_child(cell, n_bodies[i]);

    if (sc1 == sc2) {
        add_to_cell(cell->children[sc1], n_bodies, i);
    } else {
        cell->children[sc2]->index = i;  
    }
}

/* Generates the octtree for the entire system of particles */
Cell* generate_octtree(int N, Body* n_bodies) {
    // Initialize root of octtree
    Cell* root_cell = create_cell(X_BOUND, Y_BOUND, Z_BOUND);
    root_cell->index = 0;
    // cout << root_cell->n_children << endl;

    for (int i = 1; i < N; ++i) {
        Cell* cell = root_cell;

        // Find which node to add the body to
        while (cell->n_children != 0) {
            int sc = locate_child(cell, n_bodies[i]);
            cell = cell->children[sc];
        }

        add_to_cell(cell, n_bodies, i);
    }
    return root_cell;
}

/* Deletes the octtree */
void delete_octtree(Cell* cell) {
    if (cell->n_children == 0) {
        free(cell);
        return;
    }

    for (int i = 0; i < cell->n_children; ++i) {
        delete_octtree(cell->children[i]);
    }

    free(cell);
}

/* Computes the total mass and the center of mass of the current cell */
Cell* compute_cell_properties(Cell* cell, Body* n_bodies) {
    if (cell->n_children == 0) {
        if (cell->index != -1) {
            cell->center = n_bodies[cell->index];
            return cell;
        }
    } else {      
        double tx = 0, ty = 0, tz = 0;
        for (int i = 0; i < cell->n_children; ++i) {
            Cell* child = compute_cell_properties(cell->children[i], n_bodies);
            if (child != NULL) {
                (cell->center).mass += (child->center).mass;
                tx += n_bodies[child->index].px * (child->center).mass;
                ty += n_bodies[child->index].py * (child->center).mass;
                tz += n_bodies[child->index].pz * (child->center).mass;            
            }
        }
        
        // Compute center of mass
        (cell->center).px = tx / (cell->center).mass;
        (cell->center).py = ty / (cell->center).mass;
        (cell->center).pz = tz / (cell->center).mass;

        return cell;
    }
    return NULL;
}

inline double compute_distance(Body body_i, Body body_j) {
    double px_diff, py_diff, pz_diff;
    // distance in x direction
    px_diff = body_j.px - body_i.px;
    // distance in y direction
    py_diff = body_j.py - body_i.py;
    // distance in z direction
    pz_diff = body_j.pz - body_i.pz;

    // ||p_j - p_i||
    return sqrt(pow(px_diff, 2) + pow(py_diff, 2) + pow(pz_diff, 2));
}

/* Computes the force experienced between a particle and a cell */
void compute_force_from_cell(Cell* cell, int i, Body * n_bodies, double G, Force * force) {
    double px_diff, py_diff, pz_diff, factor, euclidean_distance;
    // distance in x direction
    px_diff = (cell->center).px - n_bodies[i].px;
    // distance in y direction
    py_diff = (cell->center).py - n_bodies[i].py;
    // distance in z direction
    pz_diff = (cell->center).pz - n_bodies[i].pz;

    // ||p_j - p_i||
    euclidean_distance = sqrt(pow(px_diff, 2) + pow(py_diff, 2) + pow(pz_diff, 2));

    // G * m_i * m_j / (||p_j - p_i||)^3
    factor = G * (cell->center).mass * n_bodies[i].mass / (pow(euclidean_distance, 3) + EPSILON); // + epsilon to avoid zero division
    // f_ij = factor * (p_j - p_i)
    force->fx += px_diff * factor; // force in x direction
    force->fy += py_diff * factor; // force in y direction
    force->fz += pz_diff * factor; // force in z direction 
}

/* Computes the force between the particles in the system, 
 * using the clustering-approximation for long distant forces
 */
void compute_force_from_octtree(Cell* cell, int index, Body * n_bodies, double G, Force * force) {
    if (cell->n_children == 0) {
        if (cell->index != -1 && cell->index != index) {
            compute_force_from_cell(cell, index, n_bodies, G, force);
        }
    } else {
        double d = compute_distance(n_bodies[index], cell->center);
        
        if (THETA > (cell->width / d)){ 
            // Use approximation
            compute_force_from_cell(cell, index, n_bodies, G, force);         
        } else {
            for (int i = 0; i < cell->n_children; ++i) {
                compute_force_from_octtree(cell->children[i], index, n_bodies, G, force);
            }
        }      
    }
}

// **************************************************************************
// octtree data structure and helper methods ends
// **************************************************************************

// compute force for body i based on body j where j != i
// i: body i
inline void compute_force(int i, int N, double G, Body *n_bodies, Force * force, Cell* cell){
    // reset force
    force->fx = 0;
    force->fy = 0;
    force->fz = 0;

    compute_force_from_octtree(cell, i, n_bodies, G, force);
}

// update position and velocity of a body
inline void update_body(Body * body_next, int N, double G, double TIME_DELTA, Body body_cur, Force body_force) {
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

inline void calculate(int N, int T, double G, double TIME_DELTA, Body *n_bodies) {
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

int main(int argc, char **argv) {
    uint64_t start, end;

    // n bodies
    int N;
    // n iterations
    int T;
    // frativity const
    double G;
    // time delta
    double TIME_DELTA;
    cin >> N;
    cin >> T;
    cin >> G;
    cin >> TIME_DELTA;
    
    // cout << T << " " << G << " " << TIME_DELTA << endl;

    Body n_bodies[N];
    for (int i = 0; i < N; ++i) {
        cin >> n_bodies[i].mass;
        cin >> n_bodies[i].px;
        cin >> n_bodies[i].py;
        cin >> n_bodies[i].pz;
        cin >> n_bodies[i].vx;
        cin >> n_bodies[i].vy;
        cin >> n_bodies[i].vz;
    }

    start = GetTimeStamp();
    calculate(N, T, G, TIME_DELTA, n_bodies);
    cout << "time = " << GetTimeStamp() - start << endl;

    cout << N << endl;
    for (int i = 0; i < N; ++i) {
        cout << i << ": " << n_bodies[i] << endl;
    }

    return 0;
}

// g++ -std=c++14 -O3 -o nlogn_sequential nlogn_sequential.cpp
// ./nlogn_sequential < ../body_10.data