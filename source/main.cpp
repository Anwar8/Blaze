#include <array>
#include <iostream>
#include <memory>
#include <numeric>
#include "main.hpp"
#include "Model.hpp"
#define PI 3.14159265358979323846

int main () {

    // create mesh
    real beam_length = 5.0;
    Model model;
    std::vector<coords> end_coords = {coords(0, 0, 0), coords(beam_length, 0, 0)};
    int num_divisions = 10;
    int num_elements = num_divisions;
    int num_nodes = num_elements + 1;
    model.create_line_mesh(num_divisions, end_coords);

    std::vector<unsigned> restrained_nodes = std::vector<unsigned>(num_divisions - 1);
    std::iota(restrained_nodes.begin(), restrained_nodes.end(), 2);

    // create restraints
    NodalRestraint pinned_end;
    pinned_end.assign_dofs_restraints(std::set<int>{0, 1, 2, 3, 4}); // pinned support
    pinned_end.assign_nodes_by_id(std::set<int>{1}, model.glob_mesh);

    NodalRestraint roller_support;
    roller_support.assign_dofs_restraints(std::set<int>{1, 2, 3, 4}); // roller support
    roller_support.assign_nodes_by_id(std::set<int>{num_nodes}, model.glob_mesh);

    NodalRestraint out_of_plane_restraint;  
    out_of_plane_restraint.assign_dofs_restraints(std::set<int>{1, 3, 4});
    out_of_plane_restraint.assign_nodes_by_id(restrained_nodes, model.glob_mesh);

    model.restraints.push_back(pinned_end);
    model.restraints.push_back(roller_support);
    model.restraints.push_back(out_of_plane_restraint);

    // create loads
    real y_load = -1e3;
    // buckling load is 2.58e7 N
    real buckling_load = PI*PI * (2.06e11)*(0.0004570000) / (beam_length*beam_length);
    real x_load = -1.5*buckling_load;

    model.load_manager.create_a_nodal_load_by_id({(unsigned)num_nodes}, std::set<int>{0}, std::vector<real>{x_load}, model.glob_mesh);
    model.load_manager.create_a_nodal_load_by_id({(unsigned)num_nodes/2 + 1}, std::set<int>{2}, std::vector<real>{y_load}, model.glob_mesh);

    // create a scribe and track certain DoFs
    model.scribe.track_nodes_by_id(std::set<unsigned>{(unsigned)num_nodes}, std::set<int>{0}, model.glob_mesh); 
    model.scribe.track_nodes_by_id(std::set<unsigned>{(unsigned)num_nodes/2 + 1}, std::set<int>{2}, model.glob_mesh); 
    
    // initialise restraints and loads
    model.initialise_restraints_n_loads();
    model.glob_mesh.check_nodal_loads();

    // initialise solution parameters
    real max_LF = 1;
    int nsteps = 1000;
    real tolerance = 1e-4;
    int max_iterations = 100;
    model.initialise_solution_parameters(max_LF, nsteps, tolerance, max_iterations);
    model.solve(-1);
    model.scribe.read_all_records();
}
