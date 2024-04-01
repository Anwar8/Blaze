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
    real x_load = -1e7;
    real y_load = -1e5;
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

    glob_mesh.load_node(2, 2, y_load); // load the y translation with a load for the last node (which happens to have id = 2).
    glob_mesh.load_node(2, 0, x_load); // load the x translation with a load for the last node (which happens to have id = 2).

    glob_mesh.count_dofs();
    bool converged = false;
    real tolerance = 0.002*std::max(std::abs(x_load), std::abs(y_load));
    // const std::string convergence_criterion = "norm"; // or "max" - of out of balance.
    int max_iter = 200;
    int iter = 1;

        glob_mesh.calc_global_contributions();
        assembler.assemble_global_contributions(glob_mesh);
        solver.solve_for_U(assembler);
    // begin nonlinear iterations:
    while ((iter <= max_iter) && !(converged))
    {   
        std::cout << std::endl << "-----------------------------<Iteration " << iter << ">-----------------------------" << std::endl;
        assembler.map_U_to_nodes(glob_mesh);
        glob_mesh.print_info();
        glob_mesh.update_elements_states(); // calculates internal state of strain, stress, and nodal responses
        assembler.map_elements_f_to_R(glob_mesh);
        assembler.calculate_out_of_balance();
        converged = assembler.check_convergence(tolerance);
        solver.solve_for_deltaU(assembler);
        assembler.increment_U();
        iter++;
    }

    glob_mesh.print_elements_states(true, true, true, true);
    std::cout << std::endl << "---<Analysis complete. Final iteration = " << iter << ", and out-of-balance = " << assembler.get_G_max() << ">---" << std::endl;
}
