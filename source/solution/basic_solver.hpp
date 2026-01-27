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
#ifdef WITH_BELOS
#include "BelosPseudoBlockCGSolMgr.hpp"
#include "Teuchos_ParameterList.hpp"
#else
#include "Amesos2.hpp"
#endif
#include "tpetra_wrappers.hpp"
#endif
/**
 * @brief place-holder class for solvers.
 * 
 */
class BasicSolver {
    protected:
    #ifdef WITH_MPI 
    #ifdef WITH_BELOS
    
    Teuchos::RCP<Belos::LinearProblem<scalar_type,TpetraMultiVector,operator_type>> problem_KU_P;
    Teuchos::RCP<Belos::LinearProblem<scalar_type,TpetraMultiVector,operator_type>> problem_KdU_G;

    Teuchos::RCP<Teuchos::ParameterList> belos_solver_parameters;

    Teuchos::RCP<Belos::SolverManager<scalar_type, TpetraMultiVector,operator_type> > dU_solver;
    Teuchos::RCP<Belos::SolverManager<scalar_type, TpetraMultiVector,operator_type> > U_solver;

    #else
    Teuchos::RCP<Amesos2::Solver<TpetraCrsMatrix, TpetraMultiVector>> dU_solver;
    Teuchos::RCP<Amesos2::Solver<TpetraCrsMatrix, TpetraMultiVector>> U_solver;
    #endif
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
            
            #ifdef WITH_BELOS
                belos_operator_traits::Apply(*(assembler.K), *U_rcp, *P_rcp);
                belos_operator_traits::Apply(*(assembler.K), *dU_rcp, *G_rcp);
                // Define the "problem"

                problem_KU_P = Teuchos::make_rcp<Belos::LinearProblem<scalar_type,TpetraMultiVector,operator_type>>(assembler.K, U_rcp, P_rcp);
                problem_KdU_G = Teuchos::make_rcp<Belos::LinearProblem<scalar_type,TpetraMultiVector,operator_type>>(assembler.K, dU_rcp, G_rcp);
                
                // problem_KU_P->setOperator(assembler.K);
                // problem_KdU_G->setOperator(assembler.K);

                // problem_KU_P->setRHS(P_rcp);
                // problem_KdU_G->setRHS(G_rcp);

                // problem_KU_P->setLHS(U_rcp);
                // problem_KdU_G->setLHS(dU_rcp);

                problem_KU_P->setProblem();
                problem_KdU_G->setProblem();
                // Define the solver parameters
                belos_solver_parameters = Teuchos::make_rcp<Teuchos::ParameterList>();
                int num_dofs = assembler.get_U_length();
                belos_solver_parameters->set("Maximum Iterations", 100000);       // Maximum number of iterations allowed
                belos_solver_parameters->set("Convergence Tolerance", 2e-2);         // Relative convergence tolerance 
                belos_solver_parameters->set("Output Frequency", 10000);
                belos_solver_parameters->set("Verbosity", Belos::Errors + Belos::Warnings + Belos::TimingDetails + Belos::StatusTestDetails );

                // Define the solver
                U_solver = Teuchos::make_rcp<Belos::PseudoBlockCGSolMgr<scalar_type, TpetraMultiVector, operator_type>>(problem_KU_P, belos_solver_parameters);
                dU_solver = Teuchos::make_rcp<Belos::PseudoBlockCGSolMgr<scalar_type, TpetraMultiVector, operator_type>>(problem_KdU_G, belos_solver_parameters);
                

            #else
                U_solver = Amesos2::create<TpetraCrsMatrix,TpetraMultiVector>("klu2", assembler.K, U_rcp, P_rcp);
                dU_solver = Amesos2::create<TpetraCrsMatrix,TpetraMultiVector>("klu2", assembler.K, dU_rcp, G_rcp);
            #endif
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
            #ifdef WITH_BELOS
                problem_KU_P->setProblem();
                std::cout << "Staring solving U with Belos:" << std::endl;
                // assembler.print_distributed_maths_object("U");
                // assembler.print_distributed_maths_object("P");
                U_solver->solve();
                int num_of_iterations = U_solver->getNumIters();
                std::cout << "Solver required " << num_of_iterations << " iterations." << std::endl;
                std::cout << "Finished U with Belos" << std::endl;
                // assembler.print_distributed_maths_object("U");
            #else
                U_solver->symbolicFactorization().numericFactorization().solve();
            #endif
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
            #ifdef WITH_BELOS
                problem_KdU_G->setProblem();
                std::cout << "Staring solving dU with Belos:" << std::endl;
                // assembler.print_distributed_maths_object("dU");
                // assembler.print_distributed_maths_object("G");
                dU_solver->solve();
                int num_of_iterations = dU_solver->getNumIters();
                std::cout << "Solver required " << num_of_iterations << " iterations." << std::endl;
                std::cout << "Finished dU with Belos." << std::endl;
                // assembler.print_distributed_maths_object("dU");
            #else
                dU_solver->symbolicFactorization().numericFactorization().solve();
            #endif
            #endif
        }
        
};

#endif