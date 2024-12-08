/**
 * @file SolutionProcedure.hpp
 * @brief defines the SolutionProcedure class which allows the \ref Model class to solve using a particular solution procedure such as linear, load-control, or displacement-control.
 */

#ifndef SOLUTION_PROCEDURE_HPP
#define SOLUTION_PROCEDURE_HPP

#include "maths_defaults.hpp"
#include "basic_utilities.hpp"
#include "blaze_config.hpp"
#include "global_mesh.hpp"
#include "assembler.hpp"
#include "basic_solver.hpp"
#include "LoadManager.hpp"
#include "Scribe.hpp"
#include "TimeKeeper.hpp"

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

        TimeKeeper time_keeper;
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

            time_keeper.add_timers({"all",
                                    "U_to_nodes_mapping", 
                                    "element_state_update", 
                                    "element_global_response",
                                    "assembly",
                                    "convergence_check",
                                    "dU_calculation",
                                    "material_state_update",
                                    "result_recording"});
        }

        void solve(GlobalMesh& glob_mesh, Assembler& assembler, BasicSolver& solver, LoadManager& load_manager, Scribe& scribe, int logging_frequency)
        {
            time_keeper.start_timer("all");
            while (load_factor < max_LF)
            {
                load_factor += dLF;
                #if LF_VERBOSE
                    std::cout << std::endl << "===================================[Load step " << step << " - LF = " << load_factor << "]===================================" << std::endl;
                #endif
                load_manager.increment_loads(dLF);
                bool converged = false;
                int iter = 1;

                #if VERBOSE_SLN
                    std::cout << std::endl << "Entering glob_mesh.calc_nodal_contributions_to_P()" << std::endl;
                #endif
                time_keeper.start_timer("element_global_response");
                glob_mesh.calc_nodal_contributions_to_P(); // calculates the global stiffness and load contributions from the elements and nodes, respectively. Does not assemble them.
                time_keeper.stop_timer("element_global_response");
                
                #if VERBOSE_SLN
                    std::cout << std::endl << "Entering assembler.assemble_global_P(glob_mesh)" << std::endl;
                #endif
                time_keeper.start_timer("assembly");
                assembler.assemble_global_P(glob_mesh);
                time_keeper.stop_timer("assembly");

                // begin nonlinear iterations:
                while ((iter <= max_iter) && !(converged))
                {   
                    #if LF_VERBOSE
                    std::cout << "-----------------------------------<Started: Iteration " << iter << ">-------------------------------------" << std::endl;
                    #endif
                    // solver.solve_for_U(assembler);

                    /**
                     * For the first iteration: initial mapping where none of the nodes have any displacement \f$\boldsymbol{U} = \boldsymbol{0}\f$ , and so each element \f$\boldsymbol{d} = \boldsymbol{0}\f$
                     * For each iteration after that, \f$\boldsymbol{U} = \boldsymbol{U} +  d\boldsymbol{U}\f$, and so \f$ \boldsymbol{d}\f$ maps from \f$\boldsymbol{U}\f$
                     */
                    #if VERBOSE_SLN
                    std::cout << std::endl << "Entering assembler.map_U_to_nodes(glob_mesh)" << std::endl;
                    #endif
                    time_keeper.start_timer("U_to_nodes_mapping");
                    assembler.map_U_to_nodes(glob_mesh);
                    time_keeper.stop_timer("U_to_nodes_mapping"); 
                    if (VERBOSE)
                    {
                        glob_mesh.print_info();
                    }
                    #if VERBOSE_SLN
                    std::cout << std::endl << "Entering glob_mesh.update_elements_states();" << std::endl;
                    #endif
                    time_keeper.start_timer("element_state_update");
                    glob_mesh.update_elements_states(); // calculates internal state of strain, stress, and nodal responses. 
                    time_keeper.stop_timer("element_state_update");
                    #if VERBOSE_SLN
                    std::cout << std::endl << "Entering assembler.assemble_global_K_R(glob_mesh);" << std::endl;
                    #endif
                    time_keeper.start_timer("assembly");
                    assembler.assemble_global_K_R(glob_mesh); // assembles the stiffness and load contributions into the global stiffness matrix and load vector.                    
                    time_keeper.stop_timer("assembly");
                    
                    #if VERBOSE_SLN
                    std::cout << std::endl << "Entering assembler.calculate_out_of_balance();" << std::endl;
                    #endif
                    time_keeper.start_timer("convergence_check");
                    assembler.calculate_out_of_balance();

                    #if VERBOSE_SLN
                    std::cout << std::endl << "Entering assembler.check_convergence(tolerance);" << std::endl;
                    #endif
                    converged = assembler.check_convergence(tolerance);
                    time_keeper.stop_timer("convergence_check");
                    time_keeper.start_timer("dU_calculation");
                    if (!converged)
                    {
                        #if VERBOSE_SLN
                        std::cout << std::endl << "Entering solver.solve_for_deltaU(assembler);" << std::endl;
                        #endif
                        solver.solve_for_deltaU(assembler);
                        #if VERBOSE_SLN
                        std::cout << std::endl << "Entering assembler.increment_U();" << std::endl;
                        #endif
                        assembler.increment_U();
                    }
                    time_keeper.stop_timer("dU_calculation");
                    // WARNING: this is a debugging line that MUST be removed after problem with convergence is solved.
                    // scribe.write_to_records();
                    #if LF_VERBOSE
                    if (!converged)
                    {
                        std::cout << "G_max = " << assembler.get_G_max() << " while tolerance " << tolerance << std::endl;
                        std::cout << "---------------------------------<Iteration " << iter << " Did Not Converge>--------------------------------" << std::endl;
                    } else {
                        std::cout << "-------------------------------------<Iteration " << iter << " Converged>-----------------------------------" << std::endl;
                    }
                    #endif
                    iter++;
                }
                if (VERBOSE) 
                {
                    glob_mesh.print_elements_states(true, true, true, true);
                }              
                step++;
                time_keeper.start_timer("material_state_update");
                glob_mesh.update_element_sections_starting_states();
                time_keeper.stop_timer("material_state_update");
                time_keeper.start_timer("result_recording");
                scribe.write_to_records();
                time_keeper.stop_timer("result_recording");
                if (step%logging_frequency == 0 && logging_frequency > 0)
                {
                    scribe.read_all_records();
                }
                if ((iter >= max_iter) && !(converged))
                {
                    #if LF_VERBOSE
                    std::cout << std::endl << "---<WARNING: Analysis incomplete due to convergence errors. LF = " << load_factor << ", and out-of-balance = " << assembler.get_G_max() << ">---" << std::endl;
                    #endif
                    break;
                }
            }
            #if LF_VERBOSE
            std::cout << std::endl << "---<Analysis complete. LF = " << load_factor << ", and out-of-balance = " << assembler.get_G_max() << ">---" << std::endl;
            #endif
            time_keeper.stop_timer("all");
            time_keeper.read_timers({"U_to_nodes_mapping", 
                                    "element_state_update", 
                                    "element_global_response",
                                    "assembly",
                                    "convergence_check",
                                    "dU_calculation",
                                    "material_state_update",
                                    "result_recording",
                                    "all"}, "all");
        }
};

#endif 