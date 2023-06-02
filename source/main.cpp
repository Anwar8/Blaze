#include <array>
#include <iostream>
#include "node.hpp"
#include "beam_element.hpp"
#include "mesh/global_mesh.hpp"

int main () {
    global_coords global_sys;

    global_mesh glob_mesh; 
    glob_mesh.setup_mesh("mesh/test.msh");
    glob_mesh.print_info();


    std::array<Node, 2> in_nodes;
    in_nodes[0] = Node(0.0, 0.0, 0.0);
    in_nodes[1] = Node(3.0, 0.0, 1.0);


    beam_element my_beam(in_nodes);
    my_beam.print_info();

    my_beam.calc_K();
    std::cout << "K = " << std::endl;
    std::cout << my_beam.get_K() << std::endl;

    my_beam.calc_T(global_sys.get_unit_x());
    std::cout << "T = " << std::endl;
    std::cout << my_beam.get_T() << std::endl;
    
    mat T = my_beam.get_T();
    // T(0,0) = 0.0;
    // T(1,1) = 0.0;
    // T(3,3) = 0.0;
    // T(4,4) = 0.0;
    std::cout << "Transformed K = " << std::endl;
    std::cout << T.transpose()*my_beam.get_K()*T << std::endl;
    
    
}
