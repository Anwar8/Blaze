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
        spmat K; /**< Global stiffness matrix \f$\boldsymbol{K}\f$*/
        spmat P; /**< Nodal force vector \f$\boldsymbol{P}\f$*/
        spmat R; /**< Resistance force vector \f$\boldsymbol{R}\f$*/
        spmat G; /**< Out of balance force vector  \f$\boldsymbol{G}\f$*/
        spvec U; /**< Global Nodal displacement vector \f$\boldsymbol{U}\f$ - also known as system state vector*/
    public:
        friend class BasicSolver;
        /**
         * @brief retrieves global contributions from all elements.
         * 
         * @details creates triplets by retrieving all global contributions from
         * the elements and then uses the triplets to creat the global sparse stiffness
         * matrix. Also allocates the force and displacement vectors.
         * 
         * @attention sets a force of -1e4 in one of the locations.
         * 
         * @todo add a function to read and apply forces to nodes.
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
        void map_elements_f_to_R(GlobalMesh& glob_mes);

        void calculate_out_of_balance() {
            G = K*U - R;
            std::cout << "The G (out of balance) vector is:" << std::endl << Eigen::MatrixXd(G) << std::endl;
        }
};
#endif