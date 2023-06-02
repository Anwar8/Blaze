#include <vector>
#include <tuple>
#include <string>
#include "gmsh.h"
#include "../maths_defaults.hpp"
#include "../node.hpp"

using gmsh_node_map = std::vector<std::pair<size_t, coords>>;
using gmsh_elem_map = std::vector<std::pair<size_t, std::vector<size_t>>>;

class global_mesh {
    private: 
        int nnodes = 0;
        int nelems = 0;
        int ndofs = 0;
        // std::map<size_t, std::vector<size_t>> elem_map;
        std::vector<Node> node_vector;
    public:
        void open_mesh_file(std::string const mesh_file);
        gmsh_node_map read_nodes();
        void make_nodes (gmsh_node_map node_map);
        void close_mesh_file();
        void setup_mesh(std::string const mesh_file);
        
        void print_info();


};