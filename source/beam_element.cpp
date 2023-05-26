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

void shape_function::calc_N(real x, real L) {
    N(0,0) = 1 - (x/L);
    N(0,3) = x / L;
    N(1,1) = 1 - 3*std::pow(x/L,2) + 2*std::pow(x/L,3);
    N(1,2) = x - 2*std::pow(x,2)/L + std::pow(x/L, 2)*x;
    N(1,4) = 3*std::pow(x/L, 2) - 2*std::pow(x/L, 3);
    N(1,5) = -x*(x/L) + x * std::pow(x/L,2);
}

void shape_function::calc_B(real x, real L) {
    B(0,0) = -1/L;
    B(0,3) = 1/L;
    B(1,1) = -6*std::pow(1/L,2) + 12*x*std::pow(1/L,3);
    B(1,2) = - 4/L + 6*x*std::pow(1/L, 2);
    B(1,4) = 6*std::pow(1/L, 2) - 12*x*std::pow(1/L, 3);
    B(1,5) = -2/L + 6 * x* std::pow(1/L,2);
}

void beam_element::calc_N(real x)
{
    shape_func.calc_N(x, length);
}

void beam_element::calc_B(real x)
{
    shape_func.calc_B(x, length);
}

void beam_element::calc_eps() {
    local_eps = shape_func.get_B() * local_d;
}