#ifndef NODE_HPP
#define NODE_HPP

class node {
    private:
        unsigned id = 0;
        double x, y, z, mass;
    public:
        node();
        node(double x_pos, double y_pos, double z_pos);
        void print_info();
};
#endif