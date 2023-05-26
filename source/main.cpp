#include <array>
#include <iostream>
#include "node.hpp"
#include "beam_element.hpp"

int main () {
    std::array<node, 2> in_nodes;
    in_nodes[0] = node(0.0, 0.0, 0.0);
    in_nodes[1] = node(3.0, 0.0, 0.0);

    beam_element my_beam(in_nodes);
    my_beam.print_info();
    my_beam.calc_N(1.5);
    std::cout << "N = " << std::endl << my_beam.get_N() << std::endl;
    my_beam.calc_B(1.5);
    std::cout << "B = " << std::endl << my_beam.get_B() << std::endl;
    vec d = make_xd_vec(6);
    
    std::cout << std::endl << "move 1 to the right" << std::endl;
    d(0) = 1; 
    d(3) = 1;
    my_beam.set_d(d);
    std::cout << "d = " << std::endl << my_beam.get_d() << std::endl;
    my_beam.calc_eps();
    std::cout << "eps = " << std::endl << my_beam.get_eps() << std::endl;
    
    std::cout << std::endl << "move 1 up" << std::endl;
    d.setZero();
    d(1) = 1; 
    d(4) = 1;
    my_beam.set_d(d);
    std::cout << "d = " << std::endl << my_beam.get_d() << std::endl;
    my_beam.calc_eps();
    std::cout << "eps = " << std::endl << my_beam.get_eps() << std::endl;

    std::cout << std::endl << "rotate CCW" << std::endl;       
    d.setZero();
    d(1) = -1; 
    d(2) = 2.0/3;
    d(4) = 1;
    d(5) = 2.0/3;
    my_beam.set_d(d);
    std::cout << "d = " << std::endl << my_beam.get_d() << std::endl;
    my_beam.calc_eps();
    std::cout << "eps = " << std::endl << my_beam.get_eps() << std::endl;

    std::cout << std::endl << "extend to the right by 1" << std::endl;       
    d.setZero();
    d(0) = 1; 
    my_beam.set_d(d);
    std::cout << "d = " << std::endl << my_beam.get_d() << std::endl;
    my_beam.calc_eps();
    std::cout << "eps = " << std::endl << my_beam.get_eps() << std::endl;

    std::cout << std::endl << "Bend" << std::endl;       
    d.setZero();
    d(2) = 1; 
    d(5) = -1; 
    my_beam.set_d(d);
    std::cout << "d = " << std::endl << my_beam.get_d() << std::endl;
    my_beam.calc_eps();
    std::cout << "eps = " << std::endl << my_beam.get_eps() << std::endl;

}
