#include <array>
#include <iostream>
#include <memory>
#include "node.hpp"
#include "beam_element.hpp"
#include "global_mesh.hpp"
#include "assembler.hpp"
#include "basic_solver.hpp"

int main () {
    GlobalMesh glob_mesh; 
    Assembler assembler;
    BasicSolver solver;

    glob_mesh.setup_mesh("mesh/test.msh");
    glob_mesh.count_dofs();
    int nelems = glob_mesh.get_num_elems();
    int nnodes = nelems + 1;
    glob_mesh.fix_node(1, -1);
    for (int i = 2; i <= nnodes; ++i)
    {
        glob_mesh.fix_node(i, 1);
        glob_mesh.fix_node(i, 3);
        glob_mesh.fix_node(i, 4);
    }

    glob_mesh.load_node(2, 2, -100000.); // load the y translation with a load for the last node (which happens to have id = 2).

    glob_mesh.count_dofs();
    glob_mesh.calc_global_contributions();
    assembler.assemble_global_contributions(glob_mesh);
    solver.solve_for_U(assembler);
}
