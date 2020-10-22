#ifndef __CELL_H__
#define __CELL_H__

#include "body.h"

class Cell  {
public:
    int index;                   // Index into arrays to identify particle's 
                                    // position and mass
    int n_children;              // Indicate whether cell is leaf or has 8 children
    Body center;                 // center of approximate body
                                    //  - Mass of particle of total mass of subtree
                                    //  - position of cell(cube) in space
                                    //  - position of center of mass of cell
    double width, height, depth; // Width, Height, and Depth of cell
    Cell* children[8];    // Pointers to child nodes
    // constructors
    Cell() = default;
    Cell(double width, double height, double depth)
        : width {width}, height {height}, depth {depth}, index {-1}, n_children {0}
        {
            // cout << "Cell created" << endl;
            center = Body(0.0);
        }
    
    /*
    * Generates new children for the current cell, forming a subtree. 
    * The current cell will no longer be a leaf
    */
    void generate_children() 
    {
        double w  = width / 2.0;
        double h = height / 2.0;
        double d  = depth / 2.0;

        // Cell no longer a leaf
        n_children = 8;

        // Create and initialize new children   
        for (int i = 0; i < n_children; ++i) {
            children[i] = new Cell {w, h, d};
        }

        set_location_of_children(w, h, d);
    }

    /* 
     * Locates the child to which the particle must be added 
     */
    int locate_child(const Body &body)
    {
        // Determine which child to add the body to
        if (body.px > children[6]->center.px) {
            if (body.py > children[6]->center.py) {
                if (body.pz > children[6]->center.pz) {
                    return 6;
                } else {
                    return 5;
                }
            } else {
                if (body.pz > children[6]->center.pz) {
                    return 2;
                } else {
                    return 1;
                }
            }
        } else {
            if (body.py > children[6]->center.py) {
                if (body.pz > children[6]->center.pz) {
                    return 7;
                } else {
                    return 4;
                }
            } else {
                if (body.pz > children[6]->center.pz) {
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
    void add_to_cell(Body* n_bodies, int i)
    {
        if (index == -1) {
            index = i;
            return;
        }

        generate_children();

        // The current cell's body must now be re-added to one of its children
        int sc1 = locate_child(n_bodies[index]);
        children[sc1]->index = index;

        // Locate child for new body
        int sc2 = locate_child(n_bodies[i]);

        if (sc1 == sc2) {
            children[sc1]->add_to_cell(n_bodies, i);
        } else {
            children[sc2]->index = i;
        }
    }


    /* Computes the total mass and the center of mass of the current cell */
    Cell* compute_cell_properties(Body* n_bodies) {
        if (n_children == 0) {
            if (index != -1) {
                center = n_bodies[index];
                return this;
            }
        } else {      
            double tx = 0, ty = 0, tz = 0;
            for (int i = 0; i < n_children; ++i) {
                Cell* child = children[i]->compute_cell_properties(n_bodies);

                if (child != NULL) {
                    center.mass += (child->center).mass;
                    tx += n_bodies[child->index].px * (child->center).mass;
                    ty += n_bodies[child->index].py * (child->center).mass;
                    tz += n_bodies[child->index].pz * (child->center).mass;            
                }
            }
            
            // Compute center of mass
            center.px = tx / center.mass;
            center.py = ty / center.mass;
            center.pz = tz / center.mass;

            return this;
        }
        return nullptr;
    }

private:

    void set_location_of_children(double w, double h, double d)
    {
        ((children[0])->center).set_coordinates(center);
        ((children[1])->center).set_coordinates(center, w, 0, 0);
        ((children[2])->center).set_coordinates(center, w, 0, d);
        ((children[3])->center).set_coordinates(center, 0, 0, d);
        ((children[4])->center).set_coordinates(center, 0, h, 0);
        ((children[5])->center).set_coordinates(center, w, h, 0);
        ((children[6])->center).set_coordinates(center, w, h, d);
        ((children[7])->center).set_coordinates(center, 0, h, d);
    }
};

/* Generates the octtree for the entire system of particles */
Cell* generate_octtree(int N, Body* n_bodies) {
    // Initialize root of octtree
    Cell* root_cell = new Cell {X_BOUND, Y_BOUND, Z_BOUND};
    root_cell->index = 0;
    // cout << root_cell->n_children << endl;

    for (int i = 1; i < N; ++i) {
        Cell* cell = root_cell;

        // Find which node to add the body to
        while (cell->n_children != 0) {
            int sc = cell->locate_child(n_bodies[i]);
            cell = cell->children[sc];
        }

        cell->add_to_cell(n_bodies, i);
    }
    return root_cell;
}

/* Deletes the octtree */
void delete_octtree(Cell* cell) {
    if (cell->n_children == 0) {
        delete cell;
        return;
    }

    for (int i = 0; i < cell->n_children; ++i) {
        delete_octtree(cell->children[i]);
    }

    delete cell;
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
    factor = G * (cell->center).mass * n_bodies[i].mass / (pow(euclidean_distance, 3) + EPSILON);
    
    // f_ij = factor * (p_j - p_i)
    force->fx += px_diff * factor; // force in x direction
    force->fy += py_diff * factor; // force in y direction
    force->fz += pz_diff * factor; // force in z direction 
}

/* 
 * Computes the force between the particles in the system, 
 * using the clustering-approximation for long distant forces
 */
void compute_force_from_octtree(Cell* cell, int index, Body * n_bodies, double G, Force * force) {
    if (cell->n_children == 0) {
        if (cell->index != -1 && cell->index != index) {
            compute_force_from_cell(cell, index, n_bodies, G, force);
        }
    } else {
        // double d = compute_distance(n_bodies[index], cell->center);
        double d = n_bodies[index].compute_distance(cell->center);
        
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

#endif // __CELL_H__