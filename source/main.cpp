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
    int nelems = glob_mesh.get_num_elems();
    glob_mesh.fix_node(1, -1);
    for (int i = 2; i <= nelems; ++i)
    {
        glob_mesh.fix_node(i, 1);
        glob_mesh.fix_node(i, 3);
        glob_mesh.fix_node(i, 4);
    }
    glob_mesh.count_dofs();
    glob_mesh.calc_global_contributions();
    glob_mesh.assemble_global_contributions();
    glob_mesh.solve_for_U();
}
