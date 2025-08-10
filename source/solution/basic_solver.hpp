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
#ifdef WITH_MPI
#include "Amesos2.hpp"
#include "tpetra_wrappers.hpp"
#endif
/**
 * @brief place-holder class for solvers.
 * 
 */
class BasicSolver {
    protected:
    #ifdef WITH_MPI
    // Teuchos::RCP<Amesos2::Solver<Tpetra::CrsMatrix<>, Tpetra::MultiVector<>>> dU_solver; 
    // Teuchos::RCP<Amesos2::Solver<Tpetra::CrsMatrix<>, Tpetra::MultiVector<>>> U_solver; 

    Teuchos::RCP<Amesos2::Solver<TpetraCrsMatrix, TpetraMultiVector>> dU_solver;
    Teuchos::RCP<Amesos2::Solver<TpetraCrsMatrix, TpetraMultiVector>> U_solver;

    Teuchos::RCP<TpetraMultiVector> U_rcp;
    Teuchos::RCP<TpetraMultiVector> P_rcp;
    Teuchos::RCP<TpetraMultiVector> dU_rcp;
    Teuchos::RCP<TpetraMultiVector> G_rcp;
    #else
    Eigen::SparseLU<spmat> solver;
    #endif
    public:
        /**
         * @brief creates a Teuchos::RCP that points to the memory in the Assembler that holds \f$\boldsymbol{U}\f$, \f$d\boldsymbol{U}\f$ and \f$\boldsymbol{P}\f$.
         * 
         * @param assembler the \ref Assembler object used in the \ref Model
         */
        void initialise_solver(Assembler& assembler)
        {
            #ifdef WITH_MPI
            U_rcp = Teuchos::rcpFromRef(assembler.U);
            P_rcp  = Teuchos::rcpFromRef(assembler.P);
            dU_rcp = Teuchos::rcpFromRef(assembler.dU);
            G_rcp  = Teuchos::rcpFromRef(assembler.G);
            
            // using LO = std::remove_reference_t<decltype(*assembler.K)>::local_ordinal_type;
            // using GO = std::remove_reference_t<decltype(*assembler.K)>::global_ordinal_type;
            // using ST = std::remove_reference_t<decltype(*assembler.K)>::scalar_type;

            U_solver = Amesos2::create<TpetraCrsMatrix,TpetraMultiVector>("klu2", assembler.K, U_rcp, P_rcp);
            dU_solver = Amesos2::create<TpetraCrsMatrix,TpetraMultiVector>("klu2", assembler.K, dU_rcp, G_rcp);
            #endif
        }

        /**
         * @brief solves for U using the global matrices contained in \ref Assembler; uses Eigen's SparseLU solver.
         * 
         * @param assembler 
         */
        void solve_for_U(Assembler& assembler)
        {
            #ifndef WITH_MPI
            
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
            U_solver->symbolicFactorization().numericFactorization().solve();
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
            dU_solver->symbolicFactorization().numericFactorization().solve();
            #endif
        }
        
};

#endif