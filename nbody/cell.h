#ifndef __CELL_H__
#define __CELL_H__

#include "partical.h"

class Cell  {
public:
    int index;                   // Index into arrays to identify particle's 
                                    // position and mass
    int n_children;              // Indicate whether cell is leaf or has 8 children
    Partical center;                 // center of approximate partical
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
            center = Partical(0.0);
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
    int locate_child(const Partical &partical)
    {
        // Determine which child to add the partical to
        if (partical.px > children[6]->center.px) {
            if (partical.py > children[6]->center.py) {
                if (partical.pz > children[6]->center.pz) {
                    return 6;
                } else {
                    return 5;
                }
            } else {
                if (partical.pz > children[6]->center.pz) {
                    return 2;
                } else {
                    return 1;
                }
            }
        } else {
            if (partical.py > children[6]->center.py) {
                if (partical.pz > children[6]->center.pz) {
                    return 7;
                } else {
                    return 4;
                }
            } else {
                if (partical.pz > children[6]->center.pz) {
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
    void add_to_cell(Partical* particals, int i)
    {
        if (index == -1) {
            index = i;
            return;
        }

        generate_children();

        // The current cell's partical must now be re-added to one of its children
        int sc1 = locate_child(particals[index]);
        children[sc1]->index = index;

        // Locate child for new partical
        int sc2 = locate_child(particals[i]);

        if (sc1 == sc2) {
            children[sc1]->add_to_cell(particals, i);
        } else {
            children[sc2]->index = i;
        }
    }


    /* Computes the total mass and the center of mass of the current cell */
    Cell* compute_cell_properties(Partical *particals) {
        if (n_children == 0) {
            if (index != -1) {
                center = particals[index];
                return this;
            }
        } else {      
            double tx = 0, ty = 0, tz = 0;
            for (int i = 0; i < n_children; ++i) {
                Cell* child = children[i]->compute_cell_properties(particals);

                if (child != NULL) {
                    center.mass += (child->center).mass;
                    tx += particals[child->index].px * (child->center).mass;
                    ty += particals[child->index].py * (child->center).mass;
                    tz += particals[child->index].pz * (child->center).mass;            
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

    void compute_force_from_cell(const Partical &partical, Force * force) const
    {
        double px_diff, py_diff, pz_diff, factor, d;
        // distance in x direction
        px_diff = center.px - partical.px;
        // distance in y direction
        py_diff = center.py - partical.py;
        // distance in z direction
        pz_diff = center.pz - partical.pz;

        d = partical.compute_distance(center);

        // G * m_i * m_j / (||p_j - p_i||)^3
        factor = G * center.mass * partical.mass / (pow(d, 3) + EPSILON); // + epsilon to avoid zero division
        // f_ij = factor * (p_j - p_i)
        force->fx += px_diff * factor; // force in x direction
        force->fy += py_diff * factor; // force in y direction
        force->fz += pz_diff * factor; // force in z direction 
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
Cell* generate_octtree(int N, Partical* particals) {
    // Initialize root of octtree
    Cell* root_cell = new Cell {X_BOUND, Y_BOUND, Z_BOUND};
    root_cell->index = 0;
    // cout << root_cell->n_children << endl;

    for (int i = 1; i < N; ++i) {
        Cell* cell = root_cell;

        // Find which node to add the partical to
        while (cell->n_children != 0) {
            int sc = cell->locate_child(particals[i]);
            cell = cell->children[sc];
        }

        cell->add_to_cell(particals, i);
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

/* 
 * Computes the force between the particles in the system, 
 * using the clustering-approximation for long distant forces
 */
void compute_force_from_octtree(Cell* cell, int index, Partical * particals, double G, Force * force) {
    if (cell->n_children == 0) {
        if (cell->index != -1 && cell->index != index) {
            cell->compute_force_from_cell(particals[index], force);
        }
    } else {
        // double d = compute_distance(particals[index], cell->center);
        double d = particals[index].compute_distance(cell->center);
        
        if (THETA > (cell->width / d)){ 
            // Use approximation
            cell->compute_force_from_cell(particals[index], force);
        } else {
            for (int i = 0; i < cell->n_children; ++i) {
                compute_force_from_octtree(cell->children[i], index, particals, G, force);
            }
        }      
    }
}

#endif // __CELL_H__