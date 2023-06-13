 #include <iostream>
 #include "node.hpp"
 
 
 Node::Node() : coordinates(0.0 , 0.0, 0.0), mass(0.0) {}
 Node::Node(real x_pos, real y_pos, real z_pos) : coordinates(x_pos, y_pos, z_pos), mass(0.0) {}
 
 Node::Node(int i, coords xyz) : id(i), coordinates(xyz), mass(0.0) {}

 void Node::print_info() {
    std::cout << "Node " << id << ": xyz = (" << coordinates[0] << ", " << coordinates[1] << ", " << coordinates[2] <<  "), and mass = " << mass << std::endl;
    std::cout << "There are " << std::size(connected_elements) << " connected elements. They are: ";
    print_container<std::set<int>>(connected_elements);
 }
 coords const Node::get_coords() const {
   return coordinates;
 }

void Node::fix_dof(int dof) {
  if (valid_dof(dof))
  {
    inactive_dofs.insert(dof);
  } else {
    std::cout << "ERROR: Cannot fix DoF " << dof << ". Only DoFs 0 through 5 allowed." << std::endl;
    std::exit(1);
  }
}

void Node::free_dof(int dof) {
  if (valid_dof(dof))
  {
    inactive_dofs.erase(dof);
  } else {
    std::cout << "ERROR: Cannot free DoF " << dof << ". Only DoFs 0 through 5 allowed." << std::endl;
    std::exit(1);
  }
}
