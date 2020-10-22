#include "cell.h"

void Cell::generate_children()
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

void Cell::set_location_of_children(double w, double h, double d)
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