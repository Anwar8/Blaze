#include <iostream>
#include <string>
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
    node_map.reserve(nodeTags.size());

    auto itr = coord_vec.begin();
    for (auto tag : nodeTags)
    {
        node_map.push_back(std::make_pair(tag, coords(*itr, *(itr+1), *(itr+2))));
        itr += 3;
    }
    return node_map;
}

gmsh_elem_map global_mesh::read_elements()
{
    gmsh_elem_map elem_map;

    std::vector<int> element_types;
    gmsh::model::mesh::getElementTypes(element_types);
    for (auto elem_type: element_types)
    {
        std::string type_name;
        int num_nodes_per_elem;
        // throw-away variables
        int dim, order, num_primary_nodes;
        std::vector<double> local_nodal_coords;
        // get the information about the element type. in particular:
        // (1) number of nodes, and  (2) type name
        gmsh::model::mesh::getElementProperties(elem_type, type_name, dim, 
                                                order, num_nodes_per_elem, local_nodal_coords, 
                                                num_primary_nodes);
        
        // getting the elements and nodes for each element type
        // this is because each type has the same number of nodes
        std::cout << "Mapping element type: " << type_name << std::endl;
        std::vector<std::size_t> elem_tags, node_tags;
        gmsh::model::mesh::getElementsByType(elem_type, elem_tags, node_tags);

        std::vector<std::size_t> element_nodes;
        element_nodes.reserve(num_nodes_per_elem);
        auto node_itr = node_tags.begin() ;
        for (auto elem_tag: elem_tags)
        {
            element_nodes.clear();
            for (int i = 0; i < num_nodes_per_elem; ++i)
            {
                element_nodes.push_back(*(node_itr + i));
                node_itr += i;
            }
            ++node_itr;
            elem_map.push_back(std::make_pair(elem_tag, element_nodes));
            std::cout << "added element: " << elem_tag << " with nodes: ";
            print_container<std::vector<std::size_t>>(element_nodes);
        }
    }
    return elem_map;
}
void global_mesh::make_elements (gmsh_elem_map elem_map) {
    std::vector<std::shared_ptr<Node>> elem_nodes;
    elem_nodes.reserve(2);
    for (auto element_data : elem_map)
    {
        elem_nodes.clear();
        for (auto node_id: element_data.second)
        {
            auto node = get_id_iterator<std::vector<std::shared_ptr<Node>>::iterator, std::vector<std::shared_ptr<Node>>>(node_id, node_vector);
            elem_nodes.push_back(*node);
        }
        // Basic2DBeamElement my_beam(element_data.first, elem_nodes);
        // my_beam.print_info();
        elem_vector.push_back(std::make_shared<Basic2DBeamElement>(element_data.first, elem_nodes));
    }   
}
void global_mesh::make_nodes (gmsh_node_map node_map) {
    for (auto node_data : node_map)
    {
        node_vector.push_back(std::make_shared<Node>(node_data.first, node_data.second));
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
    gmsh_elem_map elem_map = read_elements();
    nnodes = node_map.size();
    nelems = elem_map.size();
    node_vector.clear();
    node_vector.reserve(nnodes);
    elem_vector.clear();
    elem_vector.reserve(nelems);
    make_nodes(node_map);
    make_elements(elem_map);
    close_mesh_file();
}

void global_mesh::print_info()
{
    std::cout << "Mesh contains " << nelems << " elements and " << nnodes << " nodes." << std::endl;
    for (auto node: node_vector)
    {
        node->print_info();
    }
    for (auto elem: elem_vector)
    {
        elem->print_info();
    }
}

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

