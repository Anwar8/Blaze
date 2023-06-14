#include <iostream>
#include <string>
#include "beam_element.hpp"
#include "node.hpp"

Basic2DBeamElement::Basic2DBeamElement() {
}


Basic2DBeamElement::Basic2DBeamElement(std::shared_ptr<Node>& node_1, std::shared_ptr<Node>& node_2) {
    nodes[0] = node_1;
    nodes[1] = node_2;
    for (auto node : nodes) {
        node->add_connected_element(id);
    }
    calc_length();
}
Basic2DBeamElement::Basic2DBeamElement(int given_id, std::shared_ptr<Node>& node_1, std::shared_ptr<Node>& node_2) {
    id = given_id;
    nodes[0] = node_1;
    nodes[1] = node_2;
    for (auto node : nodes) {
        node->add_connected_element(id);
    }
    calc_length();
}

Basic2DBeamElement::Basic2DBeamElement(int given_id, std::vector<std::shared_ptr<Node>>& in_nodes) {
    if (std::size(in_nodes) != 2)
    {
        std::cout << "Incorrect number of element passed to create element " << id << std::endl;
        std::cout << "Received " << std::size(in_nodes) << " but expected " << 2 << std::endl; 
        std::exit(1);
    }
    id = given_id;
    nodes[0] = in_nodes[0];
    nodes[1] = in_nodes[1];
    for (auto node : in_nodes) {
        
        node->add_connected_element(id);
    }
    calc_length();
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

void BasicShapeFunction::calc_N(real x, real L) {
    N(0,0) = 1 - (x/L);
    N(0,3) = x / L;
    N(1,1) = 1 - 3*std::pow(x/L,2) + 2*std::pow(x/L,3);
    N(1,2) = x - 2*std::pow(x,2)/L + std::pow(x/L, 2)*x;
    N(1,4) = 3*std::pow(x/L, 2) - 2*std::pow(x/L, 3);
    N(1,5) = -x*(x/L) + x * std::pow(x/L,2);
}

void BasicShapeFunction::calc_B(real x, real L) {
    B(0,0) = -1/L;
    B(0,3) = 1/L;
    B(1,1) = -6*std::pow(1/L,2) + 12*x*std::pow(1/L,3);
    B(1,2) = - 4/L + 6*x*std::pow(1/L, 2);
    B(1,4) = 6*std::pow(1/L, 2) - 12*x*std::pow(1/L, 3);
    B(1,5) = -2/L + 6 * x* std::pow(1/L,2);
}

void BasicShapeFunction::calc_k(real L, BasicSection& sec) {
    std::cout << "calculating k" << std::endl;
    real A = sec.get_A();
    real E = sec.get_E();
    real I = sec.get_I();
    // Row 1
    k(0,0) = E*A/L;
    k(0,3) = -E*A/L;
    // Row 2
    k(1,1) = 12*E*I/std::pow(L,3);
    k(1,2) = 6*E*I/std::pow(L,2);
    k(1,4) = -12*E*I/std::pow(L,3);
    k(1,5) = 6*E*I/std::pow(L,2);
    // Row 3
    k(2,1) = 6*E*I/std::pow(L,2);
    k(2,2) = 4*E*I/L;
    k(2,4) = -6*E*I/std::pow(L,2);
    k(2,5) = 2*E*I/L;
    // Row 4
    k(3,0) = -E*A/L;
    k(3,3) = E*A/L;
    // Row 5
    k(4,1) = -12*E*I/std::pow(L,3);
    k(4,2) = -6*E*I/std::pow(L,2);
    k(4,4) = 12*E*I/std::pow(L,3);
    k(4,5) = -6*E*I/std::pow(L,2);
    // Row 6
    k(5,1) = 6*E*I/std::pow(L,2);
    k(5,2) = 2*E*I/L;
    k(5,4) = -6*E*I/std::pow(L,2);
    k(5,5) = 4*E*I/L;
}

void Basic2DBeamElement::calc_T(coords origin_x) {
    std::cout << "calculating T" << std::endl;
    orient.evaluate(nodes, origin_x);
}

void Basic2DBeamElement::calc_N(real x)
{
    shape_func.calc_N(x, length);
}

void Basic2DBeamElement::calc_B(real x)
{
    shape_func.calc_B(x, length);
}

void Basic2DBeamElement::calc_k()
{
    shape_func.calc_k(length, section);
}

void Basic2DBeamElement::calc_eps() {
    local_eps = shape_func.get_B() * local_d;
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

// void Basic2DBeamElement::calc_K_global() 
// {
//     calc_T();
//     calc_k();
//     mat k = orient.get_T().transpose() * shape_func.get_k() * orient.get_T();
//     K_global.clear();
//     // we have the same number of contribution as stiffness components 
//     // assuming all are non-zero!
//     K_global.reserve(k.rows() * k.cols());
//     std::vector<int> dof_map = shape_func.get_dof_map();
//     // we will first get the contribution of each node
//     int node_i = 0;
//     for (auto node: nodes)
//     {
//         std::set<int> inactive_dofs = node->get_inactive_dofs();
//         // std::set<int> active_dofs = node->get_active_dofs();
//         int node_ndof = node->get_ndof();
//         for (auto dof: active_dofs)
//         {
//             dof += node_ndof*node_i;
//         }
//         int node_id = node->get_id();
//         std::cout << "Element " << id << ", node " << node_id << std::endl;

//         for (auto node_dof: active_dofs)
//         {
//         // iterate over each element DoF
//         int elem_dof_i = 0;
//             //
//         //          MUST CHECK dof_map against inactive dofs from the node!!
//             //
            
//             for (auto elem_dof: dof_map)
//             {
//                 // check if the dof is inactive. if it is, then do nothing, else
//                 // the contributions of its element
//                 if (inactive_dofs.count(elem_dof) == 0)
//                 {
//                 // access the values by rows associated with each node
//                 // real val = k(node_i, elem_dof_i);
//                 int global_row = (node_id - 1) + node_dof;
//                 int global_column = (node_id - 1) + node_dof;
//                 ++elem_dof_i;
//                 // if this DoF is not inactive
                
//                 std::cout << "Added k(" << node_i << ", " << elem_dof_i << ")";
//                 std::cout << " to (" << global_row << ", " << global_column << ")" << std::endl;
//                 // K_global.push_back(spnz(val, global_row, global_column));

//                 }
//             }
//         }
//         // contruct the spnz from: value, row, column
//         ++node_i;
//     }

// }
void Basic2DBeamElement::calc_K_global() 
{
    std::cout << "Element " << id << " calculating its global contributions." << std::endl;
    calc_T();
    calc_k();
    std::cout << "calculating transformed k" << std::endl;
    mat k = orient.get_T().transpose() * shape_func.get_k() * orient.get_T();
    std::cout << "clearing K_global" << std::endl;
    K_global.clear();
    // we have the same number of contribution as stiffness components 
    // assuming all are non-zero!
    std::cout << "reserving space for K_global" << std::endl;
    K_global.reserve(k.rows() * k.cols());
    std::cout << "going into creating the dof map..." << std::endl;
    create_dof_map();
    
    // std::vector<int> dof_map = shape_func.get_dof_map();
    std::vector<int> dof_map = global_dof_map;
    // dof_map.insert(dof_map.end(), dof_map.begin(), dof_map.end());
    // we will first get the contribution of each node
    std::vector<int> force_in_i = dof_map;
    std::vector<int> disp_in_j = dof_map;
    int count = 1;
    int node_i_index = 0;

    for (auto node_i: nodes) {
        std::set<int> node_i_active_dofs = node_i->get_active_dofs();
        int node_i_id = node_i->get_id();
        int node_i_nz_i = node_i -> get_nz_i();
        
        int node_j_index = 0;
        for (auto node_j: nodes) {
            std::set<int> node_j_active_dofs = node_j->get_active_dofs();
            int node_j_id = node_j->get_id();
            int node_j_nz_i = node_j -> get_nz_i();
            std::cout << "Element " << id << ", nodes " << node_i_id << ", " << node_j_id << std::endl;
            for (int i = node_i_index*ndofs; i < (node_i_index + 1)*ndofs; ++i)
            {
                for (int j = node_j_index*ndofs; j < (node_j_index + 1)*ndofs; ++j)
                {
                    // -------------------------------------------------------
                    // remove this after making sure assembly is correct
                    std::cout << "i, j = " << i << ", " << j << std::endl;
                    std::cout << "node " << node_i_id << " active dofs = ";
                    print_container(node_i_active_dofs);
                    std::cout << "node " << node_j_id << " active dofs = ";
                    print_container(node_j_active_dofs);
                    // -------------------------------------------------------
                    if (node_i_active_dofs.count(force_in_i[i]) != 0 && 
                        node_j_active_dofs.count(disp_in_j[j]) != 0)
                        {
                    auto val = k(i,j);
                    int global_row = force_in_i[i] + node_i_nz_i;
                    int global_col = disp_in_j[j] + node_j_nz_i;
                    std::cout << count << ". ";
                    std::cout << "Added k(" << i << ", " << j << ")";
                    std::cout << " to (" << global_row << ", " << global_col << ")" << std::endl;
                    K_global.push_back(spnz(global_row, global_col, val));
                    }
                    ++count;   
                }
            }
        ++node_j_index;
        }
    ++node_i_index;
    }
    std::cout << std::endl << std::endl;
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

std::vector<int> map_dofs(std::vector<int> elem_dofs, std::set<int> active_dofs)
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