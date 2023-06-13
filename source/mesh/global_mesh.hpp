#include <memory>
#include <vector>
#include <tuple>
#include <string>
#include "gmsh.h"
#include "../maths_defaults.hpp"
#include "../node.hpp"
#include "../beam_element.hpp"

using gmsh_node_map = std::vector<std::pair<size_t, coords>>;
using gmsh_elem_map = std::vector<std::pair<size_t, std::vector<size_t>>>;

class global_mesh {
    private: 
        int nnodes = 0;
        int nelems = 0;
        std::vector<std::shared_ptr<Node>> node_vector;
        std::vector<std::shared_ptr<Basic2DBeamElement>> elem_vector; 

    public:
        void open_mesh_file(std::string const mesh_file);
        gmsh_node_map read_nodes();
        gmsh_elem_map read_elements();
        void make_nodes (gmsh_node_map node_map);
        void make_elements (gmsh_elem_map elem_map);
        void close_mesh_file();
        void setup_mesh(std::string const mesh_file);
        
        void print_info();
};

template <typename Iterator, typename Container>
Iterator get_id_iterator(int id, Container& a_vec)
{
    auto itr = std::begin(a_vec) + (id - 1);
    int check_id = ((*itr)->get_id());
    // this "search" is inefficient compared to other search
    // algorithms such as std::find_if or std::lower_bound
    // for the general case, but is more efficient considering
    // the average case we actually care about: a sorted
    // vector of nodes that is almost always continguous
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