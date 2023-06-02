#ifndef NODE_HPP
#define NODE_HPP
#include "maths_defaults.hpp"
class Node {
    private:
        unsigned id = 0;
        coords coordinates;
        real mass;
    public:
        Node();
        Node(real x_pos, real y_pos, real z_pos);
        Node(int i, coords xyz);
        void print_info();
        coords const get_coords() const;
        unsigned const get_id() const {return id;}
};

class global_coords {
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