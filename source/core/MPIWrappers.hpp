/**
 * @file MPIWrappers.hpp
 * @brief Header file for MPI wrapper functions.
 *
 * @details
 * This file contains the declaration of inline functions for initializing
 * and managing MPI operations. These will be mostly empty/do nothing
 * unless MPI is defined.
 */
#ifndef MPI_WRAPPERS
#define MPI_WRAPPERS

#ifndef MPI
#include "mpi.h"
inline void initialise_MPI(int& argc, char**& argv)
{
    MPI_Init(&argc, &argv);
}

inline void finalise_MPI()
{
    MPI_Finalize();
}

inline void get_my_rank(int& rank)
{
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
}

inline void get_num_ranks(int& num_ranks)
{
    MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);
}

#else

inline void initialise_MPI(int argc, char* argv[])
{

}

inline void finalise_MPI()
{

}

inline void get_my_rank(int& rank)
{
    rank = 0;
}

inline void get_num_ranks(int& num_ranks)
{
    num_ranks = 1;
}
#endif


#endif 