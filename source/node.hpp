#ifndef NODE_HPP
#define NODE_HPP
#include "maths_defaults.hpp"
class node {
    private:
        unsigned id = 0;
        real x, y, z, mass;
    public:
        node();
        node(real x_pos, real y_pos, real z_pos);
        void print_info();
};
#endif