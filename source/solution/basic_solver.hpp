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
        void solve_for_U(Assembler& assembler)
        {
            #ifndef WITH_MPI
            Eigen::SparseLU<spmat> solver;
            // Compute the ordering permutation vector from the structural pattern of A
            solver.analyzePattern(assembler.K); 
            // Compute the numerical factorization 
            solver.factorize(assembler.K); 
            //Use the factors to solve the linear system 
            
            if (solver.info() == Eigen::Success)
            {
                #if LF_VERBOSE
                std::cout << "Factorisation successful." << std::endl;
                #endif
            } else {
                std::cout << "ERROR: Factorisation unsuccessfull! Matrix is:" << std::endl;
                // convert to dense matrix to print correctly
                std::cout << Eigen::MatrixXd(assembler.K) << std::endl;

                
                std::exit(1);
            }
            assembler.U = solver.solve(assembler.P); 
            if (VERBOSE_NLB)
            {
                std::cout << "The solution is:" << std::endl << assembler.U << std::endl;
            }    
            #else
            #endif
        }
        /**
         * @brief solves for \f$\Delta \boldsymbol{U}\f$ from \f$\Delta \boldsymbol{U} = -\boldsymbol{K}^{-1} \boldsymbol{G}\f$.
         * 
         * @param assembler 
         */
        void solve_for_deltaU(Assembler& assembler)
        {
            #ifndef WITH_MPI
            Eigen::SparseLU<spmat> solver;
            // Compute the ordering permutation vector from the structural pattern of A
            solver.analyzePattern(assembler.K); 
            // Compute the numerical factorization 
            solver.factorize(assembler.K); 
            //Use the factors to solve the linear system 
            
            if (solver.info() == Eigen::Success)
            {
                #if LF_VERBOSE
                std::cout << "Factorisation successful." << std::endl;
                #endif
            } else {
                std::cout << "ERROR: Factorisation unsuccessful! Matrix is:" << std::endl;
                // convert to dense matrix to print correctly
                std::cout << Eigen::MatrixXd(assembler.K) << std::endl;

                
                std::exit(1);
            }
            assembler.dU = solver.solve(assembler.G);
            assembler.dU = -assembler.dU; 
            if (VERBOSE_NLB)
            {    
                std::cout << "dU is:" << std::endl << assembler.dU << std::endl;
            }
            #else
            #endif
        }
        
};

#endif