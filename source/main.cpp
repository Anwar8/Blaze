#include <array>
#include <iostream>
#include <memory>
#include <numeric>
#include "main.hpp"
#include "Model.hpp"
#include "ElasticPlasticMaterial.hpp"
#include "BeamColumnFiberSection.hpp"
#include "basic_utilities.hpp"

void build_an_I_section(BeamColumnFiberSection& section, ElasticPlasticMaterial& steel, real offset, real tf, real b, real tw, real h, int flange_divisions, int web_divisions)
{
    std::vector<real> areas;
    std::vector<real> ys;
    
    real y = offset;
    real area = 0.0;
    // bottom flange
    area = b*tf/(flange_divisions);
    y -= 0.5*(tf/flange_divisions);
    for (int i = 0; i < flange_divisions; ++i)
    {
        y += (tf/flange_divisions);
        ys.push_back(y);
        areas.push_back(area);
    }
    // web
    area = (h - 2*tf)*tw / web_divisions;
    y = offset + tf - 0.5*((h - 2*tf)/web_divisions);
    for (int i = flange_divisions; i < flange_divisions + web_divisions; ++i)
    {
        y += ((h - 2*tf)/web_divisions);
        ys.push_back(y);
        areas.push_back(area);
    }
    // top flange
    area = b*tf/(flange_divisions);
    y = offset + (h - tf)  - 0.5*(tf/flange_divisions);
    for (int i = flange_divisions + web_divisions; i < flange_divisions + web_divisions + flange_divisions; ++i)
    {
        y += (tf/flange_divisions);
        ys.push_back(y);
        areas.push_back(area);
    }

    section.add_fibres(&steel, areas, ys);
}

int main () {
    
    // create mesh
    real beam_length = 5.0;
    Model model;
    std::vector<coords> end_coords = {coords(0, 0, 0), coords(beam_length, 0, 0)};
    int num_divisions = 4;
    int num_elements = num_divisions;
    int num_nodes = num_elements + 1;
    
    // BasicSection sect(2.06e11, 0.0125, 0.0004570000);
    // material information
    real youngs_modulus = 2e11;
    real yield_strength = 455e6;
    real hardening_ratio = 0.01;
    ElasticPlasticMaterial steel = ElasticPlasticMaterial(youngs_modulus, yield_strength, hardening_ratio*youngs_modulus);
    // section information
    real tf = 19.6e-3;
    real tw = 11.4e-3;
    real b = 192.8e-3;
    real h = 467.2e-3;
    real moment_of_inertia = tw*pow(h - 2*tf, 3)/12 + 2*b*pow(tf,3)/12 + 2*(tf*b)*pow(0.5*h - 0.5*tf, 2); // m^4 
    BeamColumnFiberSection sect;
    build_an_I_section(sect, steel, 0.0, tf, b, tw, h, 10, 40);
    model.create_line_mesh(num_divisions, end_coords, NonlinearPlastic, sect);
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
    real moment = 1.0e4;
    real w = moment*8/(beam_length*beam_length);
    real y_load = -w*beam_length/(num_nodes - 2);
    real expected_deflection = -5 * w * std::pow(beam_length, 4)/(384*youngs_modulus*moment_of_inertia);
    // buckling load is 2.58e7 N
    // real buckling_load = PI*PI * (2.06e11)*(0.0004570000) / (beam_length*beam_length);
    // real x_load = -1.5*buckling_load;

    // model.load_manager.create_a_nodal_load_by_id({(unsigned)num_nodes}, std::set<int>{0}, std::vector<real>{x_load}, model.glob_mesh);
    model.load_manager.create_a_nodal_load_by_id({(unsigned)num_nodes/2 + 1}, std::set<int>{2}, std::vector<real>{y_load}, model.glob_mesh);

    // create a scribe and track certain DoFs
    // model.scribe.track_nodes_by_id(std::set<unsigned>{(unsigned)num_nodes}, std::set<int>{0}, model.glob_mesh); 
    model.scribe.track_nodes_by_id(std::set<unsigned>{(unsigned)num_nodes/2 + 1}, std::set<int>{2}, model.glob_mesh); 
    
    // initialise restraints and loads
    model.initialise_restraints_n_loads();
    model.glob_mesh.check_nodal_loads();

    // initialise solution parameters
    real max_LF = 1;
    int nsteps = 20;
    real tolerance = 1e-2;
    int max_iterations = 10;
    model.initialise_solution_parameters(max_LF, nsteps, tolerance, max_iterations);
    model.solve(1);
    model.scribe.read_all_records();
    std::cout << "Expected deflection is: " << expected_deflection << std::endl;
}
