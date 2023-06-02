#include <array>
#include <iostream>
#include "node.hpp"
#include "beam_element.hpp"

int main () {
    global_coords global_sys;
    std::array<node, 2> in_nodes;
    in_nodes[0] = node(0.0, 0.0, 0.0);
    in_nodes[1] = node(3.0, 0.0, 0.0);

    beam_element my_beam(in_nodes);
    my_beam.print_info();

    my_beam.calc_K();
    std::cout << "K = " << std::endl;
    std::cout << my_beam.get_K() << std::endl;

    my_beam.calc_T(global_sys.get_unit_x());
    std::cout << "T = " << std::endl;
    std::cout << my_beam.get_T() << std::endl;
    
}
