#ifndef __CELL_H__
#define __CELL_H__

#include "nlogn.h"

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
//    Cell* cell = (Cell*) malloc(sizeof(Cell));
    Cell* cell = new Cell;
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

/* 
 * Added a particle to the cell. If a particle already
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
        // free(cell);
        delete cell;
        return;
    }

    for (int i = 0; i < cell->n_children; ++i) {
        delete_octtree(cell->children[i]);
    }

    // free(cell);
    delete cell;
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
#endif // __CELL_H__