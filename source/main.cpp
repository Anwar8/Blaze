#include <array>
#include <set>
#include <vector>
#include <initializer_list>
#include <iostream>
#include <iomanip>
#include <memory>
#include <numeric>
#include "blaze_config.hpp"
#include "maths_defaults.hpp"
#include "ElementTypes.hpp"
#include "Model.hpp"
#include "ElasticPlasticMaterial.hpp"
#include "BeamColumnFiberSection.hpp"
#include "basic_utilities.hpp"
#ifdef KOKKOS
    #include <Kokkos_Core.hpp>
#endif
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

int main (int argc, char* argv[]) {
    
    // create mesh
    Model model;
    
    // BasicSection sect(2.06e11, 0.0125, 0.0004570000);
    // material information
    real youngs_modulus = 2e11;
    real yield_strength = 455e18;
    real hardening_ratio = 0.01;
    ElasticPlasticMaterial steel = ElasticPlasticMaterial(youngs_modulus, yield_strength, hardening_ratio*youngs_modulus);
    // section information
    real tf = 19.6e-3;
    real tw = 11.4e-3;
    real b = 192.8e-3;
    real h = 467.2e-3;
    real moment_of_inertia = tw*pow(h - 2*tf, 3)/12 + 2*b*pow(tf,3)/12 + 2*(tf*b)*pow(0.5*h - 0.5*tf, 2); // m^4 
    real section_area = 2*tf*b + (h - 2*tf)*tw;
    BeamColumnFiberSection sect;
    build_an_I_section(sect, steel, 0.0, tf, b, tw, h, 10, 20);
    // model.create_line_mesh(num_divisions, end_coords, NonlinearPlastic, sect);

    int nbays, nfloors, beam_divisions, column_divisions;
    real floor_height, beam_length;
    std::cout << "THERE ARE " << argc << " ARGUMENTS!!" << std::endl;
    if (argc == 3)
    {
        nbays = std::stoi(argv[1]);
        nfloors = std::stoi(argv[2]);
    }
    else
    {
        nbays = 10;
        nfloors = 5; 
    }  
    beam_divisions = 50;
    column_divisions = 35;
    floor_height = 3.5;
    beam_length = 5.0;
    model.create_frame_mesh(nbays, nfloors, beam_length, floor_height, beam_divisions, column_divisions, NonlinearPlastic, sect);
    FrameMesh the_frame = model.glob_mesh.get_frame();
    // the_frame.read_frame_size();
    std::pair<int,int> frame_size = the_frame.get_frame_size();
    std::cout << "SECTION:Model_Size" << std::endl;
    std::cout << "nbays,nfloors,beam_divisions,column_divisions,floor_height,beam_length,num_nodes,num_elements" << std::endl;
    std::cout << nbays << "," 
    << nfloors << "," << beam_divisions << "," << column_divisions << "," << floor_height << "," << beam_length << "," << frame_size.first << "," << frame_size.second << std::endl;
    std::cout << "END_SECTION:Model_Size" << std::endl;
    std::cout << "SECTION:Parallelism" << std::endl;
    read_parallelism_information();
    std::cout << "END_SECTION:Parallelism" << std::endl;

    NodalRestraint column_bases;
    column_bases.assign_dofs_restraints(std::set<int>{0, 1, 2, 3, 4, 5}); // fixed support
    column_bases.assign_nodes_by_id(the_frame.get_column_bases(), model.glob_mesh);

    NodalRestraint out_of_plane_restraint;  
    out_of_plane_restraint.assign_dofs_restraints(std::set<int>{1, 3, 4});
    out_of_plane_restraint.assign_nodes_by_id(the_frame.get_out_of_plane_nodes(), model.glob_mesh);

    model.restraints.push_back(column_bases);
    model.restraints.push_back(out_of_plane_restraint);

    std::set<size_t> loaded_nodes = the_frame.get_beam_line_node_ids(nfloors, true); 
    std::vector<unsigned> loaded_nodes_v = std::vector<unsigned>(loaded_nodes.begin(), loaded_nodes.end());
    // std::cout << "loaded nodes are: " << std::endl;    
    // print_container(loaded_nodes_v);
    model.load_manager.create_a_nodal_load_by_id(loaded_nodes_v, std::set<int>{2}, std::vector<real>{-1000}, model.glob_mesh);

    // std::set<size_t> extra_loaded_nodes = the_frame.get_beam_node_ids(12, 2, false);
    // std::vector<unsigned> extra_loaded_nodes_v = std::vector<unsigned>(extra_loaded_nodes.begin(), extra_loaded_nodes.end());
    // model.load_manager.create_a_nodal_load_by_id(extra_loaded_nodes_v, std::set<int>{2}, std::vector<real>{-500000}, model.glob_mesh);


    model.initialise_restraints_n_loads();
    model.glob_mesh.check_nodal_loads();

    // // initialise solution parameters
    real max_LF = 1;
    int nsteps = 100;
    real tolerance = 1e-2;
    int max_iterations = 10;
    model.initialise_solution_parameters(max_LF, nsteps, tolerance, max_iterations);
    #ifdef KOKKOS
        Kokkos::initialize(argc, argv);
    #endif
    model.solve(1);
    std::cout << "SECTION:Timing" << std::endl;
    model.log_timers({"U_to_nodes_mapping", 
                    "element_state_update", 
                    "element_global_response",
                    "assembly",
                    "convergence_check",
                    "dU_calculation",
                    "material_state_update",
                    "result_recording",
                    "all"});
    std::cout << "END_SECTION:Timing" << std::endl;
    #ifdef KOKKOS
        Kokkos::finalize();
    #endif
    // // model.scribe.read_all_records();
    // auto recorded_data = model.scribe.get_record_id_iterator((unsigned)num_nodes)->get_recorded_data()[2];
    
    // std::cout << std::setprecision(10); 
    // std::cout << "Computed deflection is: " << recorded_data.back() << std::endl;
    // std::cout << "Expected deflection is: " << expected_deflection << std::endl;
}
