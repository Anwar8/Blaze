#include <iostream>
#include <string>
#include "beam_element.hpp"

Basic2DBeamElement::Basic2DBeamElement() {
}

Basic2DBeamElement::Basic2DBeamElement(std::array<Node, 2> input_nodes) {
    nodes[0] = input_nodes[0];
    nodes[1] = input_nodes[1];
    calc_length();
}

void Basic2DBeamElement::print_info() {
    std::cout << elem_type << " with " << dofs << " dofs, and " << nnodes << " nodes." << std::endl;
    for (Node node_i: nodes) {
        node_i.print_info();
    }
    std::cout << "it is also of length " << length << std::endl;
}

void Basic2DBeamElement::calc_length() {
    length = (nodes[0].get_coords() - nodes[1].get_coords()).norm();
}

void BasicShapeFunction::calc_N(real x, real L) {
    N(0,0) = 1 - (x/L);
    N(0,3) = x / L;
    N(1,1) = 1 - 3*std::pow(x/L,2) + 2*std::pow(x/L,3);
    N(1,2) = x - 2*std::pow(x,2)/L + std::pow(x/L, 2)*x;
    N(1,4) = 3*std::pow(x/L, 2) - 2*std::pow(x/L, 3);
    N(1,5) = -x*(x/L) + x * std::pow(x/L,2);
}

void BasicShapeFunction::calc_B(real x, real L) {
    B(0,0) = -1/L;
    B(0,3) = 1/L;
    B(1,1) = -6*std::pow(1/L,2) + 12*x*std::pow(1/L,3);
    B(1,2) = - 4/L + 6*x*std::pow(1/L, 2);
    B(1,4) = 6*std::pow(1/L, 2) - 12*x*std::pow(1/L, 3);
    B(1,5) = -2/L + 6 * x* std::pow(1/L,2);
}

void BasicShapeFunction::calc_K(real L, BasicSection& sec) {
    real A = sec.get_A();
    real E = sec.get_E();
    real I = sec.get_I();
    // Row 1
    K(0,0) = E*A/L;
    K(0,3) = -E*A/L;
    // Row 2
    K(1,1) = 12*E*I/std::pow(L,3);
    K(1,2) = 6*E*I/std::pow(L,2);
    K(1,4) = -12*E*I/std::pow(L,3);
    K(1,5) = 6*E*I/std::pow(L,2);
    // Row 3
    K(2,1) = 6*E*I/std::pow(L,2);
    K(2,2) = 4*E*I/L;
    K(2,4) = -6*E*I/std::pow(L,2);
    K(2,5) = 2*E*I/L;
    // Row 4
    K(3,0) = -E*A/L;
    K(3,3) = E*A/L;
    // Row 5
    K(4,1) = -12*E*I/std::pow(L,3);
    K(4,2) = -6*E*I/std::pow(L,2);
    K(4,4) = 12*E*I/std::pow(L,3);
    K(4,5) = -6*E*I/std::pow(L,2);
    // Row 6
    K(5,1) = 6*E*I/std::pow(L,2);
    K(5,2) = 2*E*I/L;
    K(5,4) = -6*E*I/std::pow(L,2);
    K(5,5) = 4*E*I/L;
}

void Basic2DBeamElement::calc_T(coords origin_x) {
    orient.evaluate(nodes, origin_x);
}

void Basic2DBeamElement::calc_N(real x)
{
    shape_func.calc_N(x, length);
}

void Basic2DBeamElement::calc_B(real x)
{
    shape_func.calc_B(x, length);
}

void Basic2DBeamElement::calc_K()
{
    shape_func.calc_K(length, section);
}

void Basic2DBeamElement::calc_eps() {
    local_eps = shape_func.get_B() * local_d;
}