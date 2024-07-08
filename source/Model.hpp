/**
 * @file Model.hpp
 * @brief defines the Model class which contains the global mesh, assembler, and solver objects.
 */

#ifndef MODEL_HPP
#define MODEL_HPP
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

        void initialise_restraints_n_loads()
        {
            // apply the restraints to the global mesh and thus reduce the active freedoms. 
            for (auto restraint : restraints)
            {
                restraint.apply_restraints(glob_mesh);
            }
            // initialise all the loads from the load manager
            load_manager.initialise_loads();
        }


        void initialise_solution_parameters(real max_load_factor, int num_steps, real convergence_tolerance, int max_num_of_iterations)
        {
            solution_procedure.initialise_solution_parameters(max_load_factor, num_steps, convergence_tolerance, max_num_of_iterations);
        }

        void solve(int logging_frequency = -1)
        {
            solution_procedure.solve(glob_mesh, assembler, solver, load_manager, scribe, logging_frequency);
        }

        void create_line_mesh(int divisions, std::vector<coords> end_coords)
        {
            std::pair<NodeIdCoordsPairsVector, ElemIdNodeIdPairVector> mesh_maps = glob_mesh.map_a_line_mesh(divisions, end_coords);
            glob_mesh.setup_mesh(mesh_maps.first, mesh_maps.second);
        }
};

#endif 