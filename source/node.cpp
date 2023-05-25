 #include <iostream>
 #include "node.hpp"
 
 node::node() : x(0.0), y(0.0), z(0.0), mass(0.0) {}
 node::node(double x_pos, double y_pos, double z_pos) : x(x_pos), y(y_pos), z(z_pos), mass(0.0) {}
 void node::print_info() {
    std::cout << "x = " << x << ", y = " << y << ", z = " << z << ", and mass = " << mass << std::endl;
 }