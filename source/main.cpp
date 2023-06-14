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
    glob_mesh.count_dofs();
    glob_mesh.fix_node(1, -1);
    glob_mesh.fix_node(2, 1);
    glob_mesh.fix_node(2, 3);
    glob_mesh.fix_node(2, 4);
    glob_mesh.fix_node(3, 1);
    glob_mesh.fix_node(3, 3);
    glob_mesh.fix_node(3, 4);
    glob_mesh.count_dofs();
    glob_mesh.calc_global_contributions();
    glob_mesh.assemble_global_contributions();
    glob_mesh.solve_for_U();
}
