 #include <iostream>
 #include "node.hpp"
 
 node::node() : coordinates(0.0 , 0.0, 0.0), mass(0.0) {}
 node::node(real x_pos, real y_pos, real z_pos) : coordinates(x_pos, y_pos, z_pos), mass(0.0) {}
 void node::print_info() {
    std::cout << "xyz = (" << coordinates[0] << ", " << coordinates[1] << ", " << coordinates[2] <<  "), and mass = " << mass << std::endl;
 }
 coords const node::get_coords() const {
   return coordinates;
 }