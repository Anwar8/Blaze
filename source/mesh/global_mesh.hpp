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
Iterator get_id_iterator(int id, Container& a_vec);

