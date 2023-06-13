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
    active_dofs.erase(dof);
    calc_ndof();
  } else {
    std::cout << "ERROR: Cannot fix DoF " << dof << ". Only DoFs 0 through 5 allowed." << std::endl;
    std::exit(1);
  }
}

void Node::free_dof(int dof) {
  if (valid_dof(dof))
  {
    inactive_dofs.erase(dof);
    active_dofs.insert(dof);
    calc_ndof();
  } else {
    std::cout << "ERROR: Cannot free DoF " << dof << ". Only DoFs 0 through 5 allowed." << std::endl;
    std::exit(1);
  }
}

void Node::fix_all_dofs() {
  inactive_dofs.insert({0, 1, 2, 3, 4, 5});
  active_dofs.clear();
  calc_ndof();
}
void Node::free_all_dofs()
{
  inactive_dofs.clear();
  active_dofs.insert({0, 1, 2, 3, 4, 5});
  calc_ndof();
}

void Node::print_inactive_dofs() {
    std::cout << "Node " << id << " has " << std::size(inactive_dofs) << " inactive DoFs: ";
    print_container(inactive_dofs);
}