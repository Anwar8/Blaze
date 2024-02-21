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
 */
class Assembler {
    private:
        spmat K;
        spvec P;
        // TODO: Figure out if U should be sparse or dense!
        // vec U;
        spvec U;
    public:
        friend class BasicSolver;
        /**
         * @brief retrives global contributions from all elements.
         * 
         * @details creates triplets by retrieving all global contributions from
         * the elements and then uses the triplets to creat the global sparse stiffness
         * matrix. Also allocates the force and displacement vectors.
         * 
         * @attention sets a force of -1e4 in one of the locations.
         * 
         * @todo add a function to read and apply forces to nodes.
         * 
         * @param glob_mesh takes the global_mesh object as input to get the counters and containers for nodes and elements.
         */
        void assemble_global_contributions(GlobalMesh& glob_mesh);
};
#endif