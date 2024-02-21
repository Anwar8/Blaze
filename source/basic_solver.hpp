/**
 * @file basic_solver.hpp
 * @brief file to contain solver object and functionality to manipulating global systems of equations
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