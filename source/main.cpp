#include <array>
#include <iostream>
#include <memory>
#include "main.hpp"
#include "Model.hpp"

int main () {

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
    out_of_plane_restraint.assign_nodes_by_id(std::set<int>{2, 3, 4, 5, 6, 7, 8, 9, 10, 11}, model.glob_mesh);
    
    model.restraints.push_back(end_restraint);
    model.restraints.push_back(out_of_plane_restraint);

    // create loads
    real y_load = -1e5;
    model.load_manager.create_a_nodal_load_by_id({(unsigned)num_nodes}, std::set<int>{2}, std::vector<real>{y_load}, model.glob_mesh);

    // create a scribe and track certain DoFs
    model.scribe.track_nodes_by_id(std::set<unsigned>{(unsigned)num_nodes}, std::set<int>{2}, model.glob_mesh); 
    
    // initialise restraints and loads
    model.initialise_restraints_n_loads();
    model.glob_mesh.check_nodal_loads();

    // initialise solution parameters
    real max_LF = 1;
    int nsteps = 1;
    real tolerance = 1e-5;
    int max_iterations = 200;
    model.initialise_solution_parameters(max_LF, nsteps, tolerance, max_iterations);
    model.solve(1);
    model.scribe.read_all_records();
}
