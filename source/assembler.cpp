#include "assembler.hpp"
#include "basic_utilities.hpp"


void Assembler::assemble_global_contributions(GlobalMesh& glob_mesh) 
{
    std::vector<spnz> K_global_triplets;
    std::vector<spnz> P_global_triplets;
    std::vector<spnz> K_global_elem_triplet_contribution;
    std::vector<spnz> P_global_node_triplet_contribution;
    // TODO: this reservation is not accurate...fix it so it more closely
    // matches the number of contributions that are expected
    K_global_triplets.reserve(glob_mesh.nelems*glob_mesh.ndofs);
    P_global_triplets.reserve(glob_mesh.ndofs);

    for (auto& elem: glob_mesh.elem_vector)
    {   
        K_global_elem_triplet_contribution = elem->get_global_stiffness_triplets();
        K_global_triplets.insert(K_global_triplets.end(), K_global_elem_triplet_contribution.begin(), K_global_elem_triplet_contribution.end());
    }
    for (auto& node: glob_mesh.node_vector)
    {
        P_global_node_triplet_contribution = node->get_load_triplets();
        if (VERBOSE)
        {
        std::cout << "Assembler: Node " << node->get_id() << " with coords: "  << node->get_coords()[0] << "," << node->get_coords()[1] << "," << node->get_coords()[2] << ". its triplets are: " << std::endl;
        for (auto& triplet: P_global_node_triplet_contribution)
        {
            std::cout << "row, col, val: " << triplet.row() << "," << triplet.col() << "," << triplet.value() << std::endl;
        }
        }
        P_global_triplets.insert(P_global_triplets.end(), P_global_node_triplet_contribution.begin(), P_global_node_triplet_contribution.end());
    }
    if (VERBOSE)
    {
    std::cout << "Assembler: all triplets are: " << std::endl;
    for (auto& triplet: P_global_triplets)
    {
        std::cout << "row, col, val: " << triplet.row() << "," << triplet.col() << "," << triplet.value() << std::endl;
    }
    }
    if (VERBOSE)
    {
        std::cout << "There are " << std::size(P_global_triplets) << " P_global contributions to add up." << std::endl;
        std::cout << "There are " << std::size(K_global_triplets) << " global_stiffness_triplets contributions to add up." << std::endl;
        std::cout << "The K_global_triplets is of size " << glob_mesh.ndofs << "x" << glob_mesh.ndofs << std::endl;
    }
    K.setFromTriplets(K_global_triplets.begin(), K_global_triplets.end());
    K.makeCompressed();
    
    P.setFromTriplets(P_global_triplets.begin(), P_global_triplets.end());
    P.makeCompressed();
    if (VERBOSE_NLB)
    {
        std::cout << "The P vector is:" << std::endl << Eigen::MatrixXd(P) << std::endl;
    }
    
}

void Assembler::map_U_to_nodes(GlobalMesh& glob_mesh) 
{
    for (auto& node: glob_mesh.node_vector)
    {
        int nzi = node->get_nz_i(); // where the node displacements start in the U vector.
        int num_node_dofs = node->get_ndof(); // how many there are to loop over.
        std::set<int> node_active_dofs = node->get_active_dofs();
        int i = 0;
        
        for (auto& dof: node_active_dofs) {
            node->set_nodal_displacement(dof, U.coeff(i + nzi,0));
            ++i;
        }
    }
}

void Assembler::map_elements_f_to_R(GlobalMesh& glob_mesh) 
{
    std::vector<spnz> R_global_triplets;
    std::vector<spnz> R_global_triplets_elem_contributions;
    for (auto& elem: glob_mesh.elem_vector)
    {
        R_global_triplets_elem_contributions = elem->get_global_resistance_force_triplets();
        R_global_triplets.insert(R_global_triplets.end(), R_global_triplets_elem_contributions.begin(), R_global_triplets_elem_contributions.end());
    }
    R = make_spd_mat(glob_mesh.ndofs, 1);
    R.setFromTriplets(R_global_triplets.begin(), R_global_triplets.end());
    R.makeCompressed();
    if (VERBOSE_NLB)
    {
        std::cout << "The R vector is:" << std::endl << Eigen::MatrixXd(R) << std::endl;
        std::cout << "KU, however, is:" << std::endl << Eigen::MatrixXd(K*U) << std::endl;
        std::cout << "and P is:" << std::endl << Eigen::MatrixXd(P) << std::endl;
    }
}

bool Assembler::check_convergence(real tolerance)
{
    G_max = std::sqrt(calc_l2_norm(G));
    return G_max < tolerance;
}