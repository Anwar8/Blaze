#ifndef NODE_HPP
#define NODE_HPP
#include <set>
#include "maths_defaults.hpp"

class Node {
    private:
        unsigned id = 0;
        coords coordinates;
        real mass;
        int ndof = 6;
        std::set<int> connected_elements;
    public:
        Node();
        Node(real x_pos, real y_pos, real z_pos);
        Node(int i, coords xyz);
        void set_ndof(int dofs);
        void print_info();
        coords const get_coords() const;
        void add_connected_element(int element_id) {connected_elements.insert(element_id);}
        unsigned const get_id() const {return id;}
        
        void  set_z(real z) { coordinates[2] = z;}
};

class GlobalCoords {
    private:
        coords centroid = {0.0, 0.0 , 0.0};
        coords unit_x = {1.0, 0.0, 0.0};
        coords unit_y = {0.0, 1.0, 0.0};
        coords unit_z = {0.0, 0.0, 1.0};
    public: 
        coords get_centroid() {return centroid;}
        coords get_unit_x() {return unit_x;}
        coords get_unit_y() {return unit_y;}
        coords get_unit_z() {return unit_z;}
};

#endif