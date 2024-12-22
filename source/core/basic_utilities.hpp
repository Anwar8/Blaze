/**
 * @file basic_utilities.hpp
 * @brief basic functionality like container printing, etc.
 * 
 */

#ifndef BASIC_UTILITIES
#define BASIC_UTILITIES
#include <string>
#include <iostream>
#include <algorithm>

#ifdef KOKKOS
    #include <Kokkos_Core.hpp>
#endif
#ifdef OMP
    #include <omp.h>
#endif

/**
 * @defgroup Utility 
 * 
 * @brief utility functions used by other classes and do not generally fit elsewhere.
 * @{
 */

/**
 * @brief Get iterator for a node or element by searching for their id.
 * 
 * @details searches by relying on std::find_if from <algorithm>
 * 
 * @todo create test cases that searches when there is no node or element with the given id
 * in the searched container.
 * 
 * @tparam Iterator stl-compatible iterator corresponding to the stl compatible container.
 * @tparam Container any container with stl compatible interface.
 * @param id unique id of node or vector to search.
 * @param a_vec the container containing the nodes or the elements.
 * @return Iterator 
 */
template <typename Iterator, typename Container>
Iterator get_id_iterator(int id, Container& a_vec)
{
    Iterator itr = std::find_if(a_vec.begin(), a_vec.end(), 
        [id](const auto item) {
            return item->get_id() == id;
        });
    if (itr != a_vec.end())
    {
        return itr;
    } else 
    {
        std::cout << "could not find item with id " << id << " in vector." << std::endl;
        std::exit(1);
    }
}
/**
 * @brief Prints the contents of a container one by one.
 * 
 * @tparam T Type of data in container.
 * @param V STL-compatible container.
 */
template <typename T>
void print_container(T V)
{
  for (auto& v: V)
  {
    std::cout << v << " ";
  }
  std::cout << std::endl;
}
/** @} */ // end of Utility group

/**
 * @brief prints to the output stream the form of parallelism the program was built with, and the number of threads.
 * 
 */
void read_parallelism_information()
{
    #if defined(OMP) && !defined(KOKKOS)
    #pragma omp parallel
    {
    #pragma omp master
        std::cout <<  "parallelism,num_threads" << std::endl << "OMP," << omp_get_num_threads() << std::endl;
    }
    #elif !defined(OMP) && !defined(KOKKOS)
            std::cout <<  "parallelism,num_threads" << std::endl << "serial," << 1 << std::endl;
    #endif 

    #ifdef KOKKOS
        std::cout << "parallelism,num_threads" << std::endl << "Kokkos," << Kokkos::DefaultExecutionSpace().concurrency() << std::endl;
    #endif
}
#endif
