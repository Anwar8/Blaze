/**
 * @file basic_utilities.hpp
 * @brief basic functionality like container printing, etc.
 * 
 */

#ifndef BASIC_UTILITIES
#define BASIC_UTILITIES
#include <string>
#include <iostream>
/**
 * @defgroup Utility 
 * 
 * @brief utility functions used by other classes and do not generally fit elsewhere.
 * @{
 */

/**
 * @brief Get iterator for a node or element by searching for their id.
 * 
 * @details searches by checking if the id of the element corrsponds closely to its location
 * in the vector. if not, it goes either up or down to keep checking.
 * this "search" is inefficient compared to other search algorithms such as std::find_if
 * or std::lower_bound for the general case, but is more efficient considering the average case we 
 * actually care about: a sorted vector of nodes/elems that is almost always continguous.
 * 
 * @attention there might be a potential infinite loop in this search.
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
    auto itr = std::begin(a_vec) + (id - 1);
    int check_id = ((*itr)->get_id());
    if (check_id > id)
    {
        while (check_id != id && (itr > std::begin(a_vec)))
        {
            --itr;
            check_id = ((*itr) -> get_id());
        }
    } else if (check_id < id) {
        while (check_id != id && (itr < std::end(a_vec)))
        {
            ++itr;
            check_id = ((*itr) -> get_id());
        }
    }
    if (check_id == id)
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

#endif
