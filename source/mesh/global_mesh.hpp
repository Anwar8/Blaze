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
        int ndofs = 0;
        // std::map<size_t, std::vector<size_t>> elem_map;
        std::vector<Node> node_vector; // WARNING: does NOT work
        // if the vector is resized for any reason as all iterators
        // will be invalidated!!
        std::vector<Basic2DBeamElement> elem_vector; 

    public:
        void open_mesh_file(std::string const mesh_file);
        gmsh_node_map read_nodes();
        gmsh_elem_map read_elements();
        void make_nodes (gmsh_node_map node_map);
        void close_mesh_file();
        void setup_mesh(std::string const mesh_file);
        
        void print_info();
};

template <typename Iterator, typename Container>
Iterator get_id_iterator(int id, Container& a_vec);

void print_vector(std::vector<size_t> V);