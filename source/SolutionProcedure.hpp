/**
 * @file SolutionProcedure.hpp
 * @brief defines the SolutionProcedure class which allows the \ref Model class to solve using a particular solution procedure such as linear, load-control, or displacement-control.
 */

#ifndef SOLUTION_PROCEDURE_HPP
#define SOLUTION_PROCEDURE_HPP

#include "maths_defaults.hpp"
#include "basic_utilities.hpp"
#include "main.hpp"
#include "global_mesh.hpp"
#include "assembler.hpp"
#include "basic_solver.hpp"
#include "LoadManager.hpp"
#include "Scribe.hpp"

class SolutionProcedure
{
    protected:
        real load_factor;
        real max_LF;
        int nsteps;
        real dLF;
        int step;
        real tolerance;
        int max_iter;
    public:
        void initialise_solution_parameters(real max_load_factor, int num_steps, real convergence_tolerance, int max_num_of_iterations)
        {
            load_factor = 0;
            max_LF = max_load_factor;
            nsteps = num_steps;
            dLF = max_LF/nsteps;
            step = 1;
            tolerance = convergence_tolerance;
            max_iter = max_num_of_iterations;
        }

        void solve(GlobalMesh& glob_mesh, Assembler& assembler, BasicSolver& solver, LoadManager& load_manager, Scribe& scribe)
        {
            // should not hard-code which DOFs are tracked - that should be part of the `global_mesh` object.
            //------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            // std::vector<real> end_disp_history_y, end_disp_history_x;
            //------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            
            // I should be able to get an idea of what is loaded from the `global_mesh` object and not have to hard code the node and DoF.
            //------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            // Euler buckling load for this beam is 2.58096e7
            // real x_load = 0.0;
            // real y_load = -1e5;
            // real load_factor = 0;
            // Load is handled now by the load_manager and has already been set up elsewhere.
            //------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            
            // this really should not be part of the solution procedure - should be part of the `global_mesh` object.
            //------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            // std::pair<NodeIdCoordsPairsVector, ElemIdNodeIdPairVector> mesh_maps = glob_mesh.read_mesh_file("mesh/test.msh");
            // glob_mesh.setup_mesh(mesh_maps.first, mesh_maps.second);
            // glob_mesh.count_dofs();
            // Mesh should already be setup elsewhere.
            //------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            
            // This really needs to be handled by a different type of object, and this setup should be done elsewhere in an initialisation function.
            //------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            // int nelems = glob_mesh.get_num_elems();
            // int nnodes = nelems + 1;
            // glob_mesh.fix_node(1, -1);
            // for (int i = 2; i <= nnodes; ++i)
            // {
            //     glob_mesh.fix_node(i, 1);
            //     glob_mesh.fix_node(i, 3);
            //     glob_mesh.fix_node(i, 4);
            // }
            // glob_mesh.count_dofs();
            //------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            

            // loop over steps
            while (load_factor < max_LF)
            {
                load_factor += dLF;
                std::cout << std::endl << "-----------------------------<Load step " << step << " - LF = " << load_factor << ">-----------------------------" << std::endl;
                load_manager.increment_loads(dLF);
                bool converged = false;
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
                scribe.write_to_records();
                scribe.read_all_records();
                if ((iter >= max_iter) && !(converged))
                {
                    break;
                }
            }
            std::cout << std::endl << "---<Analysis complete. LF = " << load_factor << ", and out-of-balance = " << assembler.get_G_max() << ">---" << std::endl;
            scribe.read_all_records();
        }
};

#endif 