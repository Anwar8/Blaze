#ifndef NODE_HPP
#define NODE_HPP
#include "maths_defaults.hpp"
class node {
    private:
        unsigned id = 0;
        coords coordinates;
        real mass;
    public:
        node();
        node(real x_pos, real y_pos, real z_pos);
        void print_info();
        coords const get_coords();
};
#endif