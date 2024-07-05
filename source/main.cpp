#include <array>
#include <iostream>
#include <memory>
#include "main.hpp"
#include "Model.hpp"

int main () {
    /*
    GlobalMesh glob_mesh; 
    Assembler assembler;
    BasicSolver solver;
    
    std::vector<real> end_disp_history_y, end_disp_history_x;
    // Euler buckling load for this beam is 2.58096e7
    real x_load = 0.0;
    real y_load = -1e5;
    real load_factor = 0;
    std::pair<NodeIdCoordsPairsVector, ElemIdNodeIdPairVector> mesh_maps = glob_mesh.read_mesh_file("mesh/test.msh");
    glob_mesh.setup_mesh(mesh_maps.first, mesh_maps.second);
    glob_mesh.count_dofs();
    int nelems = glob_mesh.get_num_elems();
    int nnodes = nelems + 1;
    glob_mesh.fix_node(1, -1);
    for (int i = 2; i <= nnodes; ++i)
    {
        glob_mesh.fix_node(i, 1);
        glob_mesh.fix_node(i, 3);
        glob_mesh.fix_node(i, 4);
    }
    glob_mesh.count_dofs();

    */



    // create mesh
    Model model;
    std::vector<coords> end_coords = {coords(0, 0, 0), coords(3, 0, 0)};
    int num_divisions = 10;
    int num_elements = num_divisions;
    int num_nodes = num_elements + 1;
    model.create_line_mesh(num_divisions, end_coords);

    // create restraints
    NodalRestraint end_restraint;
    end_restraint.assign_dofs_restraints(std::set<int>{0, 1, 2, 3, 4, 5});
    end_restraint.assign_nodes_by_id(std::set<int>{1}, model.glob_mesh);

    NodalRestraint out_of_plane_restraint;  
    out_of_plane_restraint.assign_dofs_restraints(std::set<int>{1, 3, 4});
    out_of_plane_restraint.assign_nodes_by_id(std::set<int>{2, 3, 4, 5, 6, 7, 8, 9, 10}, model.glob_mesh);end_restraint.assign_nodes_by_id(std::set<int>{1}, model.glob_mesh);
    
    model.restraints.push_back(end_restraint);
    model.restraints.push_back(out_of_plane_restraint);

    // create loads
    real y_load = -1e5;
    model.load_manager.create_a_nodal_load_by_id({(unsigned)num_nodes}, std::set<int>{1}, std::vector<real>{y_load});
    
    // initialise restraints and loads
    model.initialise_restraints_n_loads();

    // initialise solution parameters
    real max_LF = 1;
    int nsteps = 100;
    real tolerance = 1e-5;
    int max_iterations = 10;
    model.initialise_solution_parameters(max_LF, nsteps, tolerance, max_iterations);
}
