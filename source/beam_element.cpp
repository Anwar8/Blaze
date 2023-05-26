#include <iostream>
#include <string>
#include "beam_element.hpp"

beam_element::beam_element() {
}

beam_element::beam_element(std::array<node, 2> input_nodes) {
    nodes[0] = input_nodes[0];
    nodes[1] = input_nodes[1];
    calc_length();
}

void beam_element::print_info() {
    std::cout << elem_type << " with " << dofs << " dofs, and " << nnodes << " nodes." << std::endl;
    for (node node_i: nodes) {
        node_i.print_info();
    }
    std::cout << "it is also of length " << length << std::endl;
}

void beam_element::calc_length() {
    length = (nodes[0].get_coords() - nodes[1].get_coords()).norm();
}