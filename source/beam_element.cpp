#include <iostream>
#include <set>
#include <string>
#include "beam_element.hpp"
#include "basic_section.hpp"
#include "basic_shape_function.hpp"
#include "basic_utilities.hpp"
#include "node.hpp"

Basic2DBeamElement::Basic2DBeamElement() {
}


Basic2DBeamElement::Basic2DBeamElement(std::shared_ptr<Node>& node_1, std::shared_ptr<Node>& node_2) {
    std::vector<std::shared_ptr<Node>> given_nodes = {node_1, node_2}; 
    initialise(0, given_nodes);
}

Basic2DBeamElement::Basic2DBeamElement(int given_id, std::shared_ptr<Node>& node_1, std::shared_ptr<Node>& node_2) {
    std::vector<std::shared_ptr<Node>> given_nodes = {node_1, node_2}; 
    initialise(given_id, given_nodes);
}

Basic2DBeamElement::Basic2DBeamElement(int given_id, std::vector<std::shared_ptr<Node>>& in_nodes) {
    initialise(given_id, in_nodes);
}

void Basic2DBeamElement::print_info() {
    std::cout << "elem " << id << " of type " <<elem_type << " with " << ndofs << " dofs, and " << nnodes << " nodes:" << std::endl;
    for (auto node_i: nodes) {
        node_i->print_info();
    }
    std::cout << "it is also of length " << length << std::endl;
}

void Basic2DBeamElement::calc_length() {
    length = (nodes[0]->get_coords() - nodes[1]->get_coords()).norm();
}

void Basic2DBeamElement::calc_T(real sec_offset, coords origin_x) {
    orient.evaluate(nodes, sec_offset, origin_x);
}

void Basic2DBeamElement::calc_N(real x)
{
    shape_func.calc_N(x, length);
}

void Basic2DBeamElement::calc_B(real x)
{
    shape_func.calc_B(x, length);
}

void Basic2DBeamElement::calc_local_constitutive_mat()
{
    real EA = section.get_E()*section.get_A();
    real EI = section.get_E()*section.get_I();
    
    local_constitutive_mat(0,0) = EA;
    local_constitutive_mat(1,0) = 0;
    local_constitutive_mat(0,1) = 0;
    local_constitutive_mat(1,1) = EI;
}
void Basic2DBeamElement::get_U_from_nodes() 
{
    std::array<real, 6> nodal_disp;
    int i = 0;
    // global_ele_U
    for (auto node: nodes)
    {
        nodal_disp = node->get_nodal_displacements();
        for (auto dof: nodal_disp)
        {
            global_ele_U(i) = dof;
            ++i;
        }
    }
}

int const Basic2DBeamElement::get_nth_node_id(int n) const {
    if (n > nnodes - 1 || n < 0)
    {
        std::cout << "Error: Requested invalid node " << n << " from element " << id << std::endl;
        std::cout << "Element has " << nnodes << " nodes." << std::endl;
        std::exit(1);
    }
    return nodes[n]->get_id();
}

void Basic2DBeamElement::map_stiffness()
{
    // local to global stiffness map: <<local_row, local_col, global_row, global_col>, ...>
    stiffness_map.clear();
    int stiffness_size = 0;
    for (auto node: nodes) 
    {
        stiffness_size += std::size(node->get_active_dofs());
    }
    stiffness_size *= stiffness_size;
    if (VERBOSE)
    {
        std::cout << "Reserved " << stiffness_size << " spaces for the stiffness_map" << std::endl;
    }
    stiffness_map.reserve(stiffness_size);
    int i = 0;
    for (auto node_i: nodes)
    {
        int j = 0;
        std::set<int> active_dofs_i = node_i->get_active_dofs();
        int nz_i_i = node_i->get_nz_i();
        for (auto node_j: nodes)
        {
            
            std::set<int> active_dofs_j = node_j->get_active_dofs();
            
            int nz_i_j = node_j->get_nz_i();
            int dof_i_index = 0;
            for (auto dof_i: active_dofs_i)
            {
                int dof_j_index = 0;
                for (auto dof_j: active_dofs_j)
                {
                    // std::cout << "i, j = " << i << ", " << j << " and their dofs are " << dof_i << ", " << dof_j << std::endl;
                    stiffness_map.push_back({6*i+dof_i, 6*j+dof_j, nz_i_i + dof_i_index, nz_i_j+dof_j_index});
                    ++dof_j_index;
                }
                ++dof_i_index;
            }
        ++j;
        }
    ++i;
    }
    if (VERBOSE)
    {
        std::cout << "Element " << id << " has " << std::size(stiffness_map) << " contributions, the stiffness map is:" << std::endl;
    
    for (auto submap: stiffness_map)
    {
        print_container(submap);
    }
    }
}
void Basic2DBeamElement::calc_K_global() 
{
    // calc_T();
    // calc_k();
    // mat k_glob = orient.get_T().transpose() * shape_func.get_k() * orient.get_T();
    global_stiffness_triplets.clear();
    // we have the same number of contribution as stiffness components 
    // assuming all are non-zero!
    global_stiffness_triplets.reserve(elem_global_stiffness.rows() * elem_global_stiffness.cols());
    for (auto kmap: stiffness_map)
    {
        real val = elem_global_stiffness(kmap[0], kmap[1]);
        global_stiffness_triplets.push_back(spnz(kmap[2], kmap[3], val));
    }
}

void Basic2DBeamElement::create_dof_map()
{
    std::cout << "creating dof map for element " << id << std::endl;
    global_dof_map.clear();
    global_dof_map.reserve(6);
    for (auto node: nodes)
    {
        std::set<int> active_dofs = node->get_active_dofs();
        std::vector<int> mapped_dofs = map_dofs(shape_func.get_dof_map(), active_dofs);
        global_dof_map.insert(global_dof_map.end(), mapped_dofs.begin(), mapped_dofs.end());
    }
}

std::vector<int> Basic2DBeamElement::map_dofs(std::vector<int> elem_dofs, std::set<int> active_dofs)
{
    std::vector<int> mapped_dofs;
    mapped_dofs.reserve(std::size(elem_dofs));
    for (auto dof: elem_dofs)
    {
        auto dof_itr = std::find(active_dofs.begin(), active_dofs.end(), dof);
        if (dof_itr == active_dofs.end())
        {
            mapped_dofs.push_back(-1);
        } else {
            mapped_dofs.push_back(std::distance(active_dofs.begin(), dof_itr));
        }
    }
    return mapped_dofs;
}
