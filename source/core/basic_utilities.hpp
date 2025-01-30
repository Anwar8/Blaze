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
 * @return Iterator that can be a vector.end() or the iterator. if vector.end() this means the id was not found.
 */
template <typename Iterator, typename Container>
Iterator get_id_iterator(int id, Container& a_vec)
{
    Iterator itr = std::find_if(a_vec.begin(), a_vec.end(), 
        [id](const auto item) {
            return item->get_id() == id;
        });
        return itr;
}

/**
 * @brief Checks whether a container of IDs are all contained in another container.
 * 
 * @details searches by relying on std::find_if from <algorithm>
 * 
 * @tparam Container any container with stl compatible interface.
 * @tparam Container any other container with stl compatible interface.
 * @param ids container of ids of nodes or vector to search.
 * @param a_vec the container containing the nodes or the elements.
 * @return Iterator 
 */
template <typename Container, typename Container2>
bool contains_id(Container& ids, Container2& a_vec)
{
    for (auto id : ids)
    {
        auto itr = std::find_if(a_vec.begin(), a_vec.end(), 
            [id](const auto item) {
                return item->get_id() == id;
            });
        if (itr == a_vec.end())
        {
            return false;
        } 
    }
    return true;
};
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
