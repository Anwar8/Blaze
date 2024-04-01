#include "basic_solver.hpp"


void BasicSolver::solve_for_U(Assembler& assembler) {
    Eigen::SparseLU<spmat> solver;
    // Compute the ordering permutation vector from the structural pattern of A
    solver.analyzePattern(assembler.K); 
    // Compute the numerical factorization 
    solver.factorize(assembler.K); 
    //Use the factors to solve the linear system 
    
    if (solver.info() == Eigen::Success)
    {
        std::cout << "Factorisation successfull." << std::endl;
    } else {
        std::cout << "ERROR: Factorisation unsuccessfull! Matrix is:" << std::endl;
        // convert to dense matrix to print correctly
        std::cout << Eigen::MatrixXd(assembler.K) << std::endl;

        
        std::exit(1);
    }
    assembler.U = solver.solve(assembler.P); 
    
    std::cout << "The solution is:" << std::endl << assembler.U << std::endl;
}
void BasicSolver::solve_for_deltaU(Assembler& assembler) {

    Eigen::SparseLU<spmat> solver;
    // Compute the ordering permutation vector from the structural pattern of A
    solver.analyzePattern(assembler.K); 
    // Compute the numerical factorization 
    solver.factorize(assembler.K); 
    //Use the factors to solve the linear system 
    
    if (solver.info() == Eigen::Success)
    {
        std::cout << "Factorisation successfull." << std::endl;
    } else {
        std::cout << "ERROR: Factorisation unsuccessfull! Matrix is:" << std::endl;
        // convert to dense matrix to print correctly
        std::cout << Eigen::MatrixXd(assembler.K) << std::endl;

        
        std::exit(1);
    }
    assembler.dU = solver.solve(assembler.G);
    assembler.dU = -assembler.dU; 
    std::cout << "dU is:" << std::endl << assembler.dU << std::endl;
}