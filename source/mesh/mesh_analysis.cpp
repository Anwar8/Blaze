#include <iostream>
#include <unordered_map>
#include "gmsh.h"

void print_vector(std::vector<size_t>);
int main(int argc, char** argv) {
    // build the geometry
    gmsh::initialize();
    gmsh::model::add("test_model");
    gmsh::open("test.msh");

    // Analysing the nodes
    std::vector<double> coords;
    std::vector<double> parametricCoords;
    std::vector<std::size_t> nodeTags;
    gmsh::model::mesh::getNodes(nodeTags, coords, parametricCoords);

    // calculate the number of nodes
    int numNodes = nodeTags.size();

    // output the results
    std::cout << "Number of nodes: " << numNodes << std::endl;

    for (auto& tag : nodeTags) {
        std::vector<double> node_coords;
        std::vector<double> param_node_coords;
        // temp1 and temp2: "dimension dim and tag tag of the entity on which the node is classified"
        int temp1, temp2;
        gmsh::model::mesh::getNode(tag, node_coords, param_node_coords, temp1, temp2);
        
        std::cout << "node " << tag << " has coordinates: " << node_coords[0] << ", " << node_coords[1] << ", " << node_coords[2] << std::endl;
    }

    // Analysing the elements
    std::vector<int> elem_types;
    std::vector<std::vector<size_t>> elem_tags, node_tags;
    
    gmsh::model::mesh::getElements(elem_types, elem_tags, node_tags);
    std::cout << "printing element tags vector:" << std::endl;
    int i = 0;
    for (auto& ele : elem_tags)
    {
        ++i;
        std::cout << "element category = " << i << std::endl;
        print_vector(ele);
    }
    std::cout << "printing node tags vector:" << std::endl;
    i = 0;
    for (auto& node : node_tags)
    {
        ++i;
        std::cout << "node category = " << i << std::endl;
        print_vector(node);
    }
    std::vector <size_t> element_tags;
    for (auto& elem_vec : elem_tags)
    {
        for (auto& tag : elem_vec)
        {
            element_tags.push_back(tag);
        }

    }
    std::unordered_map <size_t, std::vector<size_t>> elem_nodes;
    for (auto& tag : element_tags) 
    {
        std::vector<size_t> nodes;
        int element_type, dim, _tag;
        gmsh::model::mesh::getElement(tag, element_type, nodes, dim, _tag);
        elem_nodes[tag] = nodes;
    }

    for(auto const& pair: elem_nodes) {
        std::cout << "Element: " << pair.first << std::endl;
        std::cout << "nodes: ";
        print_vector(pair.second);
    }
    

    return 0;
}
void print_vector(std::vector<size_t> V) 
{
    for (auto& it = V.begin(); it != V.end(); ++it) {
        std::cout << *it;
        if (it != V.end() - 1)
        {
            std::cout << ", ";
        } else {
            std::cout << std::endl;
        }
    }
    
}