#ifndef __CELL_H__
#define __CELL_H__

#include "particle.h"

constexpr double THETA   = 1.0;        // Opening angle, for approximation in Barned hut algorithm

class Cell  {
public:
    int index;                   // Index into arrays to identify particle's 
                                    // position and mass
    int n_children;              // Indicate whether cell is leaf or has 8 children
    Particle center;                 // center of approximate particle
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
            center = default_particle();
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
        // #pragma omp parallel for num_threads(6)
        for (int i = 0; i < n_children; ++i) {
            children[i] = new Cell {w, h, d};
        }

        set_location_of_children(w, h, d);
    }

    /* 
     * Locates the child to which the particle must be added 
     */
    int locate_child(const Particle &particle)
    {
        // Determine which child to add the particle to
        if (particle.px > children[6]->center.px) {
            if (particle.py > children[6]->center.py) {
                if (particle.pz > children[6]->center.pz) {
                    return 6;
                } else {
                    return 5;
                }
            } else {
                if (particle.pz > children[6]->center.pz) {
                    return 2;
                } else {
                    return 1;
                }
            }
        } else {
            if (particle.py > children[6]->center.py) {
                if (particle.pz > children[6]->center.pz) {
                    return 7;
                } else {
                    return 4;
                }
            } else {
                if (particle.pz > children[6]->center.pz) {
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
    void add_to_cell(Particle* particles, int i)
    {
        if (index == -1) {
            index = i;
            return;
        }

        generate_children();

        // The current cell's particle must now be re-added to one of its children
        int sc1 = locate_child(particles[index]);
        children[sc1]->index = index;

        // Locate child for new particle
        int sc2 = locate_child(particles[i]);

        if (sc1 == sc2) {
            children[sc1]->add_to_cell(particles, i);
        } else {
            children[sc2]->index = i;
        }
    }


    /* 
     * Computes the total mass and the center of mass of the current cell 
     */
    Cell* generate_center(Particle *particles) {
        if (n_children == 0) {
            if (index != -1) {
                center = particles[index];
                return this;
            }
        } else {      
            double tx {0.0}, ty {0.0}, tz {0.0};
            
            for (int i = 0; i < n_children; ++i) {
                Cell* child = children[i]->generate_center(particles);

                if (child != NULL) {
                    center.mass += (child->center).mass;
                    tx += particles[child->index].px * (child->center).mass;
                    ty += particles[child->index].py * (child->center).mass;
                    tz += particles[child->index].pz * (child->center).mass;            
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

    void compute_force_from_cell(const Particle &particle, Force * force) const
    {
        double px_diff, py_diff, pz_diff, factor, d;
        // distance in x direction
        px_diff = center.px - particle.px;
        // distance in y direction
        py_diff = center.py - particle.py;
        // distance in z direction
        pz_diff = center.pz - particle.pz;

        d = compute_distance(particle, center);

        // even though we are aussing that no collide, but we still add a small number to
        // prevent 0 division
        factor = G * center.mass * particle.mass / (pow(d, 3) + EPSILON);
        force->fx += px_diff * factor; // force in x direction
        force->fy += py_diff * factor; // force in y direction
        force->fz += pz_diff * factor; // force in z direction 
    }

private:

    void set_location_of_children(double w, double h, double d)
    {
        set_coordinates((children[0])->center, center);
        set_coordinates((children[1])->center, center, w, 0, 0);
        set_coordinates((children[2])->center, center, w, 0, d);
        set_coordinates((children[3])->center, center, 0, 0, d);
        set_coordinates((children[4])->center, center, 0, h, 0);
        set_coordinates((children[5])->center, center, w, h, 0);
        set_coordinates((children[6])->center, center, w, h, d);
        set_coordinates((children[7])->center, center, 0, h, d);
    }
};


namespace BH_Octtree {
    /* 
    * Generates the octtree for the entire system of particles 
    */
    Cell* create_tree(int N, Particle* particles) {
        // Initialize root of octtree
        Cell* root_cell = new Cell {X_BOUND, Y_BOUND, Z_BOUND};
        root_cell->index = 0;

        for (int i = 1; i < N; ++i) {
            Cell* cell = root_cell;

            // Find which node to add the particle to
            while (cell->n_children != 0) {
                int sc = cell->locate_child(particles[i]);
                cell = cell->children[sc];
            }

            cell->add_to_cell(particles, i);
        }
        return root_cell;
    }

    /* 
    * Deletes the octtree and free the memory
    */
    void delete_tree(Cell* cell) {
        if (cell->n_children == 0) {
            delete cell;
            return;
        }

        for (int i = 0; i < cell->n_children; ++i) {
            delete_tree(cell->children[i]);
        }

        delete cell;
    }

    /* 
    * Computes the force between the particles in the system, 
    * using the clustering-approximation for long distant forces
    */
    void octtree_force(Cell* cell, int index, Particle * particles, Force * force) {
        if (cell->n_children == 0) {
            if (cell->index != -1 && cell->index != index) {
                cell->compute_force_from_cell(particles[index], force);
            }
        } else {
            // double d = compute_distance(particles[index], cell->center);
            double d = compute_distance(particles[index], cell->center);
            
            if (THETA > (cell->width / d)){ 
                // Use approximation
                cell->compute_force_from_cell(particles[index], force);
            } else {
                for (int i = 0; i < cell->n_children; ++i) {
                    octtree_force(cell->children[i], index, particles, force);
                }
            }      
        }
    }
}

#endif // __CELL_H__