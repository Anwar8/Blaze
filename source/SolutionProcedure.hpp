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

        void solve(GlobalMesh& glob_mesh, Assembler& assembler, BasicSolver& solver, LoadManager& load_manager, Scribe& scribe, int logging_frequency)
        {
            while (load_factor < max_LF)
            {
                load_factor += dLF;
                std::cout << std::endl << "-----------------------------<Load step " << step << " - LF = " << load_factor << ">-----------------------------" << std::endl;
                load_manager.increment_loads(dLF);
                bool converged = false;
                int iter = 1;
                
                
                // begin nonlinear iterations:
                while ((iter <= max_iter) && !(converged))
                {   
                    std::cout << std::endl << "-----------------------------<Iteration " << iter << ">-----------------------------" << std::endl;
                    
                    // solver.solve_for_U(assembler);

                    /**
                     * For the first iteration: initial mapping where none of the nodes have any displacement \f$\boldsymbol{U} = \boldsymbol{0}\f$ , and so each element \f$\boldsymbol{d} = \boldsymbol{0}\f$
                     * For each iteration after that, \f$\boldsymbol{U} = \boldsymbol{U} +  d\boldsymbol{U}\f$, and so \f$ \boldsymbol{d}\f$ maps from \f$\boldsymbol{U}\f$
                     */
                    assembler.map_U_to_nodes(glob_mesh); 
                    if (VERBOSE)
                    {
                        glob_mesh.print_info();
                    }

                    glob_mesh.update_elements_states(); // calculates internal state of strain, stress, and nodal responses. 
                    
                    glob_mesh.calc_global_contributions(); // calculates the global stiffness and load contributions from the elements and nodes, respectively. Does not assemble them.
                    assembler.assemble_global_contributions(glob_mesh); // assembles the stiffness and load contributions into the global stiffness matrix and load vector.

                    assembler.map_elements_f_to_R(glob_mesh);

                    assembler.calculate_out_of_balance();
                    converged = assembler.check_convergence(tolerance);
                    if (!converged)
                    {
                        solver.solve_for_deltaU(assembler);
                        assembler.increment_U();
                    }
                    // WARNING: this is a debugging line that MUST be removed after problem with convergence is solved.
                    // scribe.write_to_records();
                    std::cout << std::endl << "-----------------------------<Completed: Iteration " << iter << ">-----------------------------" << std::endl;
                    iter++;
                }
                if (VERBOSE) 
                {
                    glob_mesh.print_elements_states(true, true, true, true);
                }              
                step++;
                scribe.write_to_records();
                if (step%logging_frequency == 0 && logging_frequency > 0)
                {
                    scribe.read_all_records();
                }
                if ((iter >= max_iter) && !(converged))
                {
                    std::cout << std::endl << "---<WARNING: Analysis incomplete due to convergence errors. LF = " << load_factor << ", and out-of-balance = " << assembler.get_G_max() << ">---" << std::endl;
                    break;
                }
            }
            std::cout << std::endl << "---<Analysis complete. LF = " << load_factor << ", and out-of-balance = " << assembler.get_G_max() << ">---" << std::endl;
            if (logging_frequency > 0)
            {
                scribe.read_all_records();
            }
            
        }
};

#endif 