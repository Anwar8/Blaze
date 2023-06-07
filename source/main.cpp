#include <array>
#include <iostream>
#include <memory>
#include "node.hpp"
#include "beam_element.hpp"
#include "mesh/global_mesh.hpp"

int main () {
    GlobalCoords global_sys;

    global_mesh glob_mesh; 
    glob_mesh.setup_mesh("mesh/test.msh");
    glob_mesh.print_info();
    std::shared_ptr<Node> in_nodes_1 = std::make_shared<Node>(0.0, 0.0, 0.0);
    std::shared_ptr<Node> in_nodes_2 = std::make_shared<Node>(3.0, 0.0, 1.0);



    Basic2DBeamElement my_beam(in_nodes_1, in_nodes_2);
    my_beam.print_info();
    my_beam.move_nodes_up(66.66);
    std::cout << "moved nodes up..." << std::endl << std::endl;
    std::cout << "checking beam info..." << std::endl;
    my_beam.print_info();
    std::cout << "checking parent array..." << std::endl;
    in_nodes_1->print_info();
    in_nodes_2->print_info();


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
