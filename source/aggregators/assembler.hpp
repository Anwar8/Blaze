/**
 * @file assembler.hpp
 * @brief assembler object and its related functionality
 * 
 */
#ifndef ASSEMBLER
#define ASSEMBLER

#include "global_mesh.hpp"
#ifdef KOKKOS
    #include <Kokkos_Core.hpp>
#endif
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

        std::vector<spnz> K_global_triplets; /** The container for the triplets that are used for assembling the global stiffness matrix \f$ \boldsymbol{K}\f$ from element contributions.*/
        std::vector<spnz> R_global_triplets; /** The container for the triplets that are used for assembling the global resistance vector \f$ \boldsymbol{R}\f$ from element contributions.*/
        std::vector<spnz> P_global_triplets;  /** The container for the triplets that are used for assembling the global load matrix \f$ \boldsymbol{P}\f$ from nodal contributions.*/
    public:
        friend class BasicSolver;
        /**
         * @brief initialises the K, P, U, and dU sparse matrices and vectors to the sizes that correspond to the mesh being used.
         * 
         * @param glob_mesh  the global_mesh object which contains information about the number of nodes and degrees of freedom.
         */
        void initialise_global_matrices(GlobalMesh& glob_mesh) {
            K = make_spd_mat(glob_mesh.ndofs, glob_mesh.ndofs);
            R = make_spd_mat(glob_mesh.ndofs, 1);
            P = make_spd_mat(glob_mesh.ndofs, 1);
            U = make_spd_vec(glob_mesh.ndofs);
            dU = make_spd_vec(glob_mesh.ndofs);

            K_global_triplets.reserve(glob_mesh.nelems*36);
            R_global_triplets.reserve(glob_mesh.ndofs);
            P_global_triplets.reserve(glob_mesh.ndofs);
        }
        /**
         * @brief retrieves global load contributions from all nodes.
         * 
         * @details creates triplets by retrieving all global contributions from
         * the nodes and then uses the triplets to create the global load vector \f$ \boldsymbol{P}\f$
         * Note that as per the documentation of Eigen: "The initial contents of *this is destroyed." 
         * when calling [setFromTriplets](https://eigen.tuxfamily.org/dox-devel/classEigen_1_1SparseMatrix.html#a8f09e3597f37aa8861599260af6a53e0).
         * So this means that this function does not add to the existing force vector but recreates it from the contributions of the nodes.
         * 
         * @param glob_mesh takes the global_mesh object as input to get the counters and containers for nodes and elements.
         */
        void assemble_global_P(GlobalMesh& glob_mesh)
        {
            P_global_triplets.clear();

            for (auto& node: glob_mesh.node_vector)
            {
                node->insert_load_triplets(this->P_global_triplets);
            }
            if (VERBOSE)
            {
                std::cout << "Assembler: all triplets are: " << std::endl;
                for (auto& triplet: P_global_triplets)
                {
                    std::cout << "row, col, val: " << triplet.row() << "," << triplet.col() << "," << triplet.value() << std::endl;
                }
                std::cout << "There are " << std::size(P_global_triplets) << " P_global contributions to add up." << std::endl;
            }
            P.setFromTriplets(P_global_triplets.begin(), P_global_triplets.end());
            P.makeCompressed();
            if (VERBOSE_NLB)
            {
                std::cout << "The P vector is:" << std::endl << Eigen::MatrixXd(P) << std::endl;
            }
        }
        /**
         * @brief retrieves global contributions to stiffness and resistance from all elements. For resistance, this function maps local element nodal forces \f$\boldsymbol{f}\f$ to the resistance vector \f$\boldsymbol{R}\f$ \ref R.
         * 
         * @details creates triplets by retrieving all global contributions from
         * the elements and then uses the triplets to create the global sparse stiffness
         * matrix, and the resistance vector. Note that as per the 
         * documentation of Eigen: "The initial contents of *this is destroyed." when calling [setFromTriplets](https://eigen.tuxfamily.org/dox-devel/classEigen_1_1SparseMatrix.html#a8f09e3597f37aa8861599260af6a53e0).
         * So this means that this function does not add to the existing matrix but recreates it from the contributions of the elements.
         * @todo add a function to calculate constrained nodes reactions.
         * 
         * @param glob_mesh takes the global_mesh object as input to get the counters and containers for nodes and elements.
         */
        void assemble_global_K_R(GlobalMesh& glob_mesh)
        {
            // clear the members but keep the size reservation unchanged.
            K_global_triplets.clear();
            R_global_triplets.clear();
            for (auto& elem: glob_mesh.elem_vector)
            {   
                elem->insert_global_stiffness_triplets(K_global_triplets);
                elem->insert_global_resistance_force_triplets(R_global_triplets);
            }
            K.setFromTriplets(K_global_triplets.begin(), K_global_triplets.end());
            K.makeCompressed();
            R.setFromTriplets(R_global_triplets.begin(), R_global_triplets.end());
            R.makeCompressed();

            if (VERBOSE)
            {
                std::cout << "There are " << std::size(K_global_triplets) << " global_stiffness_triplets contributions to add up." << std::endl;
                std::cout << "The K_global_triplets is of size " << glob_mesh.ndofs << "x" << glob_mesh.ndofs << std::endl;
            }
            if (VERBOSE_NLB)
            {
                std::cout << "The R vector is:" << std::endl << Eigen::MatrixXd(R) << std::endl;
                std::cout << "KU, however, is:" << std::endl << Eigen::MatrixXd(K*U) << std::endl;
                std::cout << "and P is:" << std::endl << Eigen::MatrixXd(P) << std::endl;
            }
        }

        /**
         * @brief maps state vector U back to nodes.
         * 
         * @param glob_mesh takes the global_mesh object as input to get the counters and containers for nodes and elements.
         */
        void map_U_to_nodes(GlobalMesh& glob_mesh)
        {
            std::vector<std::shared_ptr<Node>>* nodes = &glob_mesh.node_vector;
            #ifdef KOKKOS
                Kokkos::parallel_for( "Assembler::map_U_to_nodes", glob_mesh.node_vector.size(), KOKKOS_LAMBDA (int i) {
                    int nzi = (*nodes)[i]->get_nz_i(); // where the node displacements start in the U vector.
                    int num_node_dofs = (*nodes)[i]->get_ndof(); // how many there are to loop over.
                    std::set<int> node_active_dofs = (*nodes)[i]->get_active_dofs();
                    int j = 0;
                    for (auto& dof: node_active_dofs) {
                        (*nodes)[i]->set_nodal_displacement(dof, U.coeff(j + nzi,0));
                        ++j;
                    }
                });
            #else
                #pragma omp parallel for
                for (auto& node: glob_mesh.node_vector)
                {
                    int nzi = node->get_nz_i(); // where the node displacements start in the U vector.
                    int num_node_dofs = node->get_ndof(); // how many there are to loop over.
                    std::set<int> node_active_dofs = node->get_active_dofs();
                    int i = 0;
                    
                    for (auto& dof: node_active_dofs) {
                        node->set_nodal_displacement(dof, U.coeff(i + nzi,0));
                        ++i;
                    }
                }
            #endif
        }

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
        bool check_convergence(real tolerance)
        {
            G_max = std::sqrt(calc_l2_norm(G));
            return G_max < tolerance;
        }

        realx2 get_G_max() {return G_max;}
};
#endif