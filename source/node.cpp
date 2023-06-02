 #include <iostream>
 #include "node.hpp"
 
 Node::Node() : coordinates(0.0 , 0.0, 0.0), mass(0.0) {}
 Node::Node(real x_pos, real y_pos, real z_pos) : coordinates(x_pos, y_pos, z_pos), mass(0.0) {}
 
 Node::Node(int i, coords xyz) : id(i), coordinates(xyz), mass(0.0) {}
 
 void Node::print_info() {
    std::cout << "xyz = (" << coordinates[0] << ", " << coordinates[1] << ", " << coordinates[2] <<  "), and mass = " << mass << std::endl;
 }
 coords const Node::get_coords() const {
   return coordinates;
 }