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

class Model
{
    public:
        GlobalMesh glob_mesh; 
        Assembler assembler;
        BasicSolver solver;
        SolutionProcedure solution_procedure;

};

#endif 