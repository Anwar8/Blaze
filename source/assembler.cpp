#include "assembler.hpp"


void Assembler::assemble_global_contributions(GlobalMesh& glob_mesh) 
{
    std::vector<spnz> K_global;
    std::vector<spnz> K_global_elem;
    // TODO: this reservation is not accurate...fix it so it more closely
    // matches the number of contributions that are expected
    K_global.reserve(glob_mesh.nelems*glob_mesh.ndofs);
    for (auto elem: glob_mesh.elem_vector)
    {   
        K_global_elem = elem->get_K_global();
        K_global.insert(K_global.end(), K_global_elem.begin(), K_global_elem.end());
    }
    K = make_spd_mat(glob_mesh.ndofs, glob_mesh.ndofs);
    P = make_spd_vec(glob_mesh.ndofs);
    // TODO: Decide whether U is sparse or dense
    // U = make_xd_vec(glob_mesh.ndofs);
    U = make_spd_vec(glob_mesh.ndofs);
    std::cout << "There are " << std::size(K_global) << " contributions to add up." << std::endl;
    std::cout << "The K_global is of size " << glob_mesh.ndofs << "x" << glob_mesh.ndofs << std::endl;
    K.setFromTriplets(K_global.begin(), K_global.end());
    K.makeCompressed();
    std::cout << "Setting a force of -1e4 N on node in y direction." << std::endl;
    P.insert(1) = -1e4;
}
