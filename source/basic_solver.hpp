/**
 * @file basic_solver.hpp
 * @brief file to contain solver object and functionality to manipulating global systems of equations
 * 
 * @todo need to map state vector (displacements U) back to nodes
 * @todo need to have state of calculating element strains, stresses and resistances based on nodal displacements
 * @todo state recording and I/O necessary
 * @todo wrap entire solution procedure so that solver can have an "perform_iterative_step" function. Eigen3-like or PyTorch-like solver API
 * 
 */
#ifndef BASIC_SOLVER
#define BASIC_SOLVER
#include "assembler.hpp"

/**
 * @brief place-holder class for solvers.
 * 
 */
class BasicSolver {
    public:
        /**
         * @brief solves for U using the global matrices contained in \ref Assembler; uses Eigen's SparseLU solver.
         * 
         * @param assembler 
         */
        void solve_for_U(Assembler& assembler);
};

#endif