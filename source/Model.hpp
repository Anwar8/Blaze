/**
 * @file Model.hpp
 * @brief defines the Model class which contains the global mesh, assembler, and solver objects.
 */

#ifndef MODEL_HPP
#define MODEL_HPP
#include "ElementTypes.hpp"
#include "global_mesh.hpp"
#include "assembler.hpp"
#include "basic_solver.hpp"
#include "SolutionProcedure.hpp"
#include "LoadManager.hpp"
#include "Scribe.hpp"
#include "NodalRestraint.hpp"

class Model
{
    public:
        GlobalMesh glob_mesh; 
        Assembler assembler;
        BasicSolver solver;
        SolutionProcedure solution_procedure;
        LoadManager load_manager;
        Scribe scribe;
        std::vector<NodalRestraint> restraints;

        /**
         * @brief Initialise restraints and loads for the model.
         *
         * This function applies the restraints to the global mesh, reducing the active freedoms.
         * It then initialises all the loads from the load manager. Based on the restraints,
         * the global matrices are initialised and element-stiffness mapping is established. 
         * 
         * This docstring was generated by Copilot but is correct.
         */
        void initialise_restraints_n_loads()
        {
            // apply the restraints to the global mesh and thus reduce the active freedoms. 
            #ifndef WITH_MPI
                for (auto& restraint : restraints)
                {
                    restraint.apply_restraints(glob_mesh);
                }
            #else
                for (auto& restraint : restraints)
                {
                    restraint.apply_restraints();
                }
                glob_mesh.count_and_exchange_distributed_dofs();
                glob_mesh.find_max_num_stiffness_contributions();
            #endif
            // initialise all the loads from the load manager
            load_manager.initialise_loads();
            glob_mesh.map_element_stiffnesses();

            glob_mesh.calc_nodal_contributions_to_P();


            // based on the restraints and loads, we can now initialise the global matrices and establish element-stiffness mapping.
            // start with mapping stiffnesses as these are needed for initialising the sparse matrices (the stiffness matrix).
            assembler.initialise_global_vectors(glob_mesh);
            assembler.assemble_global_P(glob_mesh);

            assembler.map_U_to_nodes(glob_mesh);
            glob_mesh.update_elements_states();
            assembler.initialise_stiffness_matrix(glob_mesh);
            // std::cout << "The initial stiffness matrix K is:" << std::endl;
            // assembler.print_distributed_maths_object("K", Teuchos::VERB_EXTREME);
            assembler.assemble_global_K_R(glob_mesh);
        }


        void initialise_solution_parameters(real max_load_factor, int num_steps, real convergence_tolerance, int max_num_of_iterations)
        {
            solution_procedure.initialise_parallel_timer(glob_mesh.get_mesh_rank(), glob_mesh.get_mesh_num_ranks());
            solution_procedure.initialise_solution_parameters(max_load_factor, num_steps, convergence_tolerance, max_num_of_iterations);
            solver.initialise_solver(assembler);
        }

        void solve(int logging_frequency = -1)
        {
            solution_procedure.solve(glob_mesh, assembler, solver, load_manager, scribe, logging_frequency);
        }



        void create_line_mesh(int divisions, std::vector<coords> end_coords, ElementType elem_type, BeamColumnFiberSection& sect)
        {
            glob_mesh.create_line_mesh(divisions, end_coords, elem_type, sect);
        }
        void create_line_mesh(int divisions, std::vector<coords> end_coords, ElementType elem_type, BasicSection& sect)
        {
            glob_mesh.create_line_mesh(divisions, end_coords, elem_type, sect);
        }
        void create_frame_mesh(int nbays, int nfloors, real bay_length, real floor_height, int beam_divisions, int column_divisions, ElementType elem_type, BeamColumnFiberSection& sect)
        {
            glob_mesh.create_frame_mesh(nbays, nfloors, bay_length, floor_height, beam_divisions, column_divisions, elem_type, sect);
        }



        void create_distributed_line_mesh(int divisions, std::vector<coords> end_coords, ElementType elem_type, BeamColumnFiberSection& sect)
        {
            glob_mesh.create_distributed_line_mesh(divisions, end_coords, elem_type, sect);
        }
        void create_distributed_line_mesh(int divisions, std::vector<coords> end_coords, ElementType elem_type, BasicSection& sect)
        {
            glob_mesh.create_distributed_line_mesh(divisions, end_coords, elem_type, sect);
        }
        void create_distributed_frame_mesh(int nbays, int nfloors, real bay_length, real floor_height, int beam_divisions, int column_divisions, ElementType elem_type, BeamColumnFiberSection& sect)
        {
            glob_mesh.create_distributed_frame_mesh(nbays, nfloors, bay_length, floor_height, beam_divisions, column_divisions, elem_type, sect);
        }


        void read_all_records()
        {
            scribe.read_all_records();
        }

        void log_timers(std::vector<std::string> timers_names)
        {
            solution_procedure.log_timers(timers_names);
        }
        
        /**
         * @brief calls `log_parallel_timers` from the \ref TimeKeeper stored in \ref SolutionProcedure.
         * @warning cannot be 
         * 
         * @param timers_names 
         */
        void log_parallel_timers(std::vector<std::string> timers_names)
        {
            solution_procedure.log_parallel_timers(timers_names);
        }

        void read_timers(std::vector<std::string> timers_names, std::string reference_timer = "")
        {
            solution_procedure.read_timers(timers_names, reference_timer);
        }
};

#endif 