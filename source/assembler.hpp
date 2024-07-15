/**
 * @file assembler.hpp
 * @brief assembler object and its related functionality
 * 
 */
#ifndef ASSEMBLER
#define ASSEMBLER

#include "global_mesh.hpp"
/**
 * @brief helps assemble the global matrices.
 * 
 * @details supposed to help with organising the code by separating 
 * assembly; right now it is not helpful and only contains one big
 * function and redundunt member variables that appeared in \ref GlobalMesh.
 * 
 * @todo improve teh naming scheme of the variables. Right now they are either very
 * vague such as "K_glob" (now became K_glob_triplets) when it is actually referring to triplet vector, or way too
 * long such as "K_global_elem_triplet_contribution"
 * 
 * @todo break down the assembler into load and stiffness assembly steps in stead of combined together in a big mess.
 * 
 */
class Assembler {
    private:
        spmat K; /**< Global stiffness matrix \f$\boldsymbol{K}\f$.*/
        spmat P; /**< Nodal force vector \f$\boldsymbol{P}\f$.*/
        spmat R; /**< Resistance force vector \f$\boldsymbol{R}\f$.*/
        spmat G; /**< Out of balance force vector  \f$\boldsymbol{G}\f$.*/
        realx2 G_max; /**< Max out-of-balance force retrieved by taking the square-root of the l2 norm of \f$ \sqrt{\norm{\boldsymbol{G}}}\f$.*/
        spvec U; /**< Global Nodal displacement vector \f$\boldsymbol{U}\f$ - also known as system state vector.*/
        spvec dU; /**< Change in global Nodal displacement vector \f$\Delta\boldsymbol{U}\f$ - also known as system state vector increment.*/
    public:
        friend class BasicSolver;
        /**
         * @brief initialises the K, P, U, and dU sparse matrices and vectors to the sizes that correspond to the mesh being used.
         * 
         * @param glob_mesh  the global_mesh object which contains information about the number of nodes and degrees of freedom.
         */
        void initialise_global_matrices(GlobalMesh& glob_mesh) {
            K = make_spd_mat(glob_mesh.ndofs, glob_mesh.ndofs);
            P = make_spd_mat(glob_mesh.ndofs, 1);
            U = make_spd_vec(glob_mesh.ndofs);
            dU = make_spd_vec(glob_mesh.ndofs);
        }
        /**
         * @brief retrieves global contributions from all elements.
         * 
         * @details creates triplets by retrieving all global contributions from
         * the elements and then uses the triplets to creat the global sparse stiffness
         * matrix. Also allocates the force and displacement vectors. Note that as per the 
         * documentation of Eigen: "The initial contents of *this is destroyed." when calling [setFromTriplets](https://eigen.tuxfamily.org/dox-devel/classEigen_1_1SparseMatrix.html#a8f09e3597f37aa8861599260af6a53e0).
         * So this means that this function does not add to the existing matrix but recreates it from the contributions of the elements.
         * @todo add a function to calculate constrained nodes reactions.
         * 
         * @param glob_mesh takes the global_mesh object as input to get the counters and containers for nodes and elements.
         */
        void assemble_global_contributions(GlobalMesh& glob_mesh);

        /**
         * @brief maps state vector U back to nodes.
         * 
         * @param glob_mesh takes the global_mesh object as input to get the counters and containers for nodes and elements.
         */
        void map_U_to_nodes(GlobalMesh& glob_mesh);

        /**
         * @brief maps local element nodal forces \f$\boldsymbol{f}\f$ to the resistance vector \f$\boldsymbol{R}\f$ \ref R.
         * 
         * @param glob_mesh takes the global_mesh object as input to get the counters and containers for nodes and elements.
         */
        void map_elements_f_to_R(GlobalMesh& glob_mesh);
        /**
         * @brief calculates out of balance forces from \f$\boldsymbol{G} =  \boldsymbol{R} - \boldsymbol{P}\f$.
         * 
         */
        void calculate_out_of_balance() {
            G = R - P;
            if (VERBOSE_NLB)
            {
                std::cout << "The G (out of balance) vector is:" << std::endl << Eigen::MatrixXd(G) << std::endl;
            }
        }
        /**
         * @brief updates \f$\boldsymbol{U}\f$ by incrementing with \f$ \Delta \boldsymbol{U}\f$.
         * 
         */
        void increment_U() {
            if (VERBOSE_NLB)
                std::cout << "U before update is " << std::endl << U << std::endl;
            U += dU;
            if (VERBOSE_NLB)
                std::cout << "U after update is " << std::endl << U << std::endl;
        }

        /**
         * @brief checks if the maximum square-root of the norm of out-of-balance is smaller than a tolerance.
         * 
         * @param tolerance maximum allowed out of balance force.
         * @return true if converged (i.e. tolerance more than max out-of-balance)
         * @return false if not converged (i.e. tolerance more than out-of-balance)
         */
        bool check_convergence(real tolerance);

        realx2 get_G_max() {return G_max;}
};
#endif