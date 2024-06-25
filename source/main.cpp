#include <array>
#include <iostream>
#include <memory>
#include "main.hpp"
#include "basic_utilities.hpp"
#include "node.hpp"
#include "beam_element.hpp"
#include "global_mesh.hpp"
#include "assembler.hpp"
#include "basic_solver.hpp"

int main () {
    GlobalMesh glob_mesh; 
    Assembler assembler;
    BasicSolver solver;
    
    std::vector<real> end_disp_history_y, end_disp_history_x;
    // Euler buckling load for this beam is 2.58096e7
    real x_load = 0.0;
    real y_load = -1e5;
    real load_factor = 0;
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
    glob_mesh.count_dofs();


    real max_LF = 1;
    int nsteps = 100;
    real dLF = max_LF/nsteps;
    // int recording_interval = 10;
    
    int step = 1;
    real LF = 0;

    // loop over steps
    while (LF < max_LF)
    {
        LF += dLF;
        std::cout << std::endl << "-----------------------------<Load step " << step << " - LF = " << LF << ">-----------------------------" << std::endl;
        if (step != 1)
        {
            glob_mesh.increment_node_load(2, 2, dLF*y_load); // load the y translation with a load for the last node (which happens to have id = 2).
            glob_mesh.increment_node_load(2, 0, dLF*x_load); // load the x translation with a load for the last node (which happens to have id = 2).
        } else {
            glob_mesh.load_node(2, 2, dLF*y_load); // load the y translation with a load for the last node (which happens to have id = 2).
            glob_mesh.load_node(2, 0, dLF*x_load); // load the x translation with a load for the last node (which happens to have id = 2).
        }
        
        bool converged = false;
        // real tolerance = 0.00002*std::max(std::abs(x_load), std::abs(y_load));
        real tolerance = 1000;
        // const std::string convergence_criterion = "norm"; // or "max" - of out of balance.
        int max_iter = 20;
        int iter = 1;

        glob_mesh.calc_global_contributions();
        assembler.assemble_global_contributions(glob_mesh);
        solver.solve_for_U(assembler);
        // begin nonlinear iterations:
        while ((iter <= max_iter) && !(converged))
        {   
            std::cout << std::endl << "-----------------------------<Iteration " << iter << ">-----------------------------" << std::endl;
            assembler.map_U_to_nodes(glob_mesh);
            if (VERBOSE)
            {
                glob_mesh.print_info();
            }
            glob_mesh.update_elements_states(); // calculates internal state of strain, stress, and nodal responses
            assembler.map_elements_f_to_R(glob_mesh);
            assembler.calculate_out_of_balance();
            converged = assembler.check_convergence(tolerance);
            solver.solve_for_deltaU(assembler);
            assembler.increment_U();
            std::cout << std::endl << "-----------------------------<Completed: Iteration " << iter << ">-----------------------------" << std::endl;
            iter++;
        }
        if (VERBOSE) 
        {
            glob_mesh.print_elements_states(true, true, true, true);
        }
        
        step++;
        
        // if (!(step%recording_interval))
        // {
            glob_mesh.track_nodal_dof(2, 2, end_disp_history_y);
            glob_mesh.track_nodal_dof(2, 0, end_disp_history_x);
        // }
        std::cout << "LF = " << LF << " and node 2 DoF 2 (U2) = " << *(end_disp_history_y.end()-1) << std::endl;
        std::cout << "LF = " << LF << " and node 2 DoF 0 (U1) = " << *(end_disp_history_x.end()-1) << std::endl;
        if ((iter >= max_iter) && !(converged))
        {
            break;
        }
    }
    std::cout << std::endl << "---<Analysis complete. LF = " << LF << ", and out-of-balance = " << assembler.get_G_max() << ">---" << std::endl;
    std::cout << "nodal displacement history for DoF 2 of node 2 is: " << std::endl;
    print_container(end_disp_history_y);
    std::cout << "nodal displacement history for DoF 0 of node 2 is: " << std::endl;
    print_container(end_disp_history_x);
    
}
