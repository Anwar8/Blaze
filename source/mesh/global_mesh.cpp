#include "global_mesh.hpp"

void global_mesh::open_mesh_file(std::string const mesh_file) {
    gmsh::initialize();
    gmsh::open(mesh_file);
}

gmsh_node_map global_mesh::read_nodes() {
    std::vector<double> coord_vec;
    std::vector<double> parametricCoords;
    std::vector<std::size_t> nodeTags;

    gmsh::model::mesh::getNodes(nodeTags, coord_vec, parametricCoords);
    
    gmsh_node_map node_map;
    auto itr = coord_vec.begin();
    for (auto tag : nodeTags)
    {
        node_map.push_back(std::make_pair(tag, coords(*itr, *(itr+1), *(itr+2))));
    }
    return node_map;
}

void global_mesh::make_nodes (gmsh_node_map node_map) {
    for (auto node_data : node_map)
    {
        node_vector.push_back(Node(node_data.first, node_data.second));
    }
}

void global_mesh::close_mesh_file()
{
    gmsh::finalize();
}

void global_mesh::setup_mesh(std::string const mesh_file) 
{
    open_mesh_file(mesh_file);
    gmsh_node_map node_map = read_nodes();
    nnodes = node_map.size();
    node_vector.clear();
    make_nodes(node_map);
    close_mesh_file();
}