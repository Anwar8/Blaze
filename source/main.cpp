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
#include "MPIWrappers.hpp"
#include "tpetra_wrappers.hpp"
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

/**
 * @brief options that user can change; sets default values.
 * 
 */
struct InputOptions {
    int nbays = 10;
    int nfloors = 5;
    int flange_divisions = 10;
    int web_divisions = 40;

    real tf = 19.6e-3;
    real tw = 11.4e-3;
    real b = 192.8e-3;
    real h = 467.2e-3;
    real yield_strength = 455e6;
    real hardening_ratio = 0.01;
    real youngs_modulus = 2e11;

    real udl = -3000; 

    int beam_divisions = 50;
    int column_divisions = 35;
    real floor_height = 3.5;
    real beam_length = 5.0;

    real max_LF = 1;
    int nsteps = 10;
    real tolerance = 1e-2;
    int max_iterations = 10;

    ElementType element_type = LinearElastic;
    BasicSection basic_sect;
    BeamColumnFiberSection fibre_sect;

    real get_nodal_load()
    {
        return udl * beam_length / (beam_divisions - 1);
    }

    void set_element_type(std::string element_type_str) {
        if (element_type_str == "NonlinearPlastic")
            element_type =  NonlinearPlastic;
        else if (element_type_str == "NonlinearElastic")
            element_type =  NonlinearElastic;
        else if (element_type_str == "LinearElastic")
            element_type = LinearElastic;
        else
        {
            std::cout << "Element string is " << element_type_str << ", while only accept LinearElastic, NonlinearElastic, and NonlinearPlastic" << std::endl;
            exit(-1);
        }
    }

    void print_summary_csv(int num_nodes, int num_elements) const {
    std::cout << "nbays,nfloors,beam_length,floor_height,beam_divisions,column_divisions,element_type,tf,tw,b,h,yield_strength,hardening_ratio,youngs_modulus,flange_divisions,web_divisions,udl,max_LF,nsteps,tolerance,max_iterations,num_nodes,num_elements" << std::endl;
    std::cout << nbays << ","
                << nfloors << ","
                << beam_length << ","
                << floor_height << ","
                << beam_divisions << ","
                << column_divisions << ","
                << element_type << ","
                << tf << ","
                << tw << ","
                << b << ","
                << h << ","
                << yield_strength << ","
                << hardening_ratio << ","
                << youngs_modulus << ","
                << flange_divisions << ","
                << web_divisions << ","
                << udl << ","
                << max_LF << ","
                << nsteps << ","
                << tolerance << ","
                << max_iterations << ","
                << num_nodes << ","
                << num_elements << std::endl;
    }
};

/**
 * @brief populates the input options for Blaze.
 * 
 * @param argc 
 * @param argv 
 * @return InputOptions 
 */
InputOptions parse_input(int argc, char* argv[]) {
    InputOptions opts;
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--nbays" && i + 1 < argc) {
            opts.nbays = std::stoi(argv[++i]);
        } else if (arg == "--nfloors" && i + 1 < argc) {
            opts.nfloors = std::stoi(argv[++i]);
        } else if (arg == "--flange_divisions" && i + 1 < argc) {
            opts.flange_divisions = std::stoi(argv[++i]);
        } else if (arg == "--web_divisions" && i + 1 < argc) {
            opts.web_divisions = std::stoi(argv[++i]);
        } else if (arg == "--elem_type" && i + 1 < argc) {
            opts.set_element_type(argv[++i]);
        } else if (arg == "--udl" && i + 1 < argc) {
            opts.udl = std::stod(argv[++i]);
        } else if (arg == "--beam_divisions" && i + 1 < argc) {
            opts.beam_divisions = std::stoi(argv[++i]);
        } else if (arg == "--column_divisions" && i + 1 < argc) {
            opts.column_divisions = std::stoi(argv[++i]);
        } else if (arg == "--floor_height" && i + 1 < argc) {
            opts.floor_height = std::stod(argv[++i]);
        } else if (arg == "--beam_length" && i + 1 < argc) {
            opts.beam_length = std::stod(argv[++i]);
        } else if (arg == "--max_LF" && i + 1 < argc) {
            opts.max_LF = std::stod(argv[++i]);
        } else if (arg == "--nsteps" && i + 1 < argc) {
            opts.nsteps = std::stoi(argv[++i]);
        } else if (arg == "--tolerance" && i + 1 < argc) {
            opts.tolerance = std::stod(argv[++i]);
        } else if (arg == "--max_iterations" && i + 1 < argc) {
            opts.max_iterations = std::stoi(argv[++i]);
        } else if (arg == "--tf" && i + 1 < argc) {
            opts.tf = std::stod(argv[++i]);
        } else if (arg == "--tw" && i + 1 < argc) {
            opts.tw = std::stod(argv[++i]);
        } else if (arg == "--b" && i + 1 < argc) {
            opts.b = std::stod(argv[++i]);
        } else if (arg == "--h" && i + 1 < argc) {
            opts.h = std::stod(argv[++i]);
        } else if (arg == "--yield_strength" && i + 1 < argc) {
            opts.yield_strength = std::stod(argv[++i]);
        } else if (arg == "--youngs_modulus" && i + 1 < argc) {
            opts.youngs_modulus = std::stod(argv[++i]);
        } else if (arg == "--hardening_ratio" && i + 1 < argc) {
            opts.hardening_ratio = std::stod(argv[++i]);
        } else {
            std::cout << "ignoring unknown option " << arg << "." << std::endl;
        }

    }
    return opts;
}
int main (int argc, char* argv[]) {

    #ifdef KOKKOS
        Kokkos::initialize(argc, argv);
    #endif
    
    #ifdef WITH_MPI
    Tpetra::ScopeGuard tpetraScope (&argc, &argv);
    #endif
    {
    //input handling
    #ifdef WITH_MPI
    Teuchos::RCP<const Teuchos::Comm<int>> comm = Teuchos::rcp(new Teuchos::MpiComm<int>(MPI_COMM_WORLD));
    #endif
    TimeKeeper time_keeper;
    Model model;
    int rank = -1;
    int num_ranks = -1;
    get_my_rank(rank);
    get_num_ranks(num_ranks);
    time_keeper.initialise_parallel_keeper(rank, num_ranks);
    std::vector<std::string> timers_names =  {"mesh_setup", "bc_load_records", "initialisation", "solution", "all"};
    time_keeper.add_timers(timers_names);
    time_keeper.start_timer("all");
    InputOptions input_options = parse_input(argc, argv);
    time_keeper.start_timer("mesh_setup");
    // material information
    ElasticPlasticMaterial steel = ElasticPlasticMaterial(input_options.youngs_modulus, input_options.yield_strength, input_options.hardening_ratio*input_options.youngs_modulus);

    // section information
    real moment_of_inertia = input_options.tw*pow(input_options.h - 2*input_options.tf, 3)/12 + 2*input_options.b*pow(input_options.tf,3)/12 + 2*(input_options.tf*input_options.b)*pow(0.5*input_options.h - 0.5*input_options.tf, 2); // m^4 
    real section_area = 2*input_options.tf*input_options.b + (input_options.h - 2*input_options.tf)*input_options.tw;
    build_an_I_section(input_options.fibre_sect, steel, 0.0, input_options.tf, input_options.b, input_options.tw, input_options.h, input_options.flange_divisions, input_options.web_divisions);
    input_options.basic_sect = BasicSection(input_options.youngs_modulus, section_area, moment_of_inertia);
    
    
    // Definition of the mesh:
    if (input_options.element_type == NonlinearPlastic)
    {
        model.create_distributed_frame_mesh(input_options.nbays, input_options.nfloors, input_options.beam_length, input_options.floor_height, input_options.beam_divisions, input_options.column_divisions, input_options.element_type, input_options.fibre_sect);
    }
    else
    {
        model.create_distributed_frame_mesh(input_options.nbays, input_options.nfloors, input_options.beam_length, input_options.floor_height, input_options.beam_divisions, input_options.column_divisions, input_options.element_type, input_options.basic_sect);
    }
    time_keeper.stop_timer("mesh_setup");
    // BCs, Loads, and Records
    time_keeper.start_timer("bc_load_records");
    FrameMesh the_frame = model.glob_mesh.get_frame();
    std::pair<int,int> frame_size = the_frame.get_frame_size();

    time_keeper.stop_timer("bc_load_records"); // we don't want the logging to take up time
    if (rank == 0)
    {
        input_options.print_summary_csv(frame_size.first, frame_size.second);
    }
    time_keeper.start_timer("bc_load_records");
    
    // Boundary conditions
    NodalRestraint column_bases;
    column_bases.assign_dofs_restraints(std::set<int>{0, 1, 2, 3, 4, 5}); // fixed support
    column_bases.assign_distributed_nodes_by_record_id(the_frame.get_column_bases(), model.glob_mesh);

    NodalRestraint out_of_plane_restraint;  
    out_of_plane_restraint.assign_dofs_restraints(std::set<int>{1, 3, 4});
    out_of_plane_restraint.assign_distributed_nodes_by_record_id(the_frame.get_out_of_plane_nodes(), model.glob_mesh);

    model.restraints.push_back(column_bases);
    model.restraints.push_back(out_of_plane_restraint);

    // Loads
    std::set<unsigned> loaded_nodes = the_frame.get_all_beam_line_node_ids(false); 
    std::vector<unsigned> loaded_nodes_v = std::vector<unsigned>(loaded_nodes.begin(), loaded_nodes.end());
    model.load_manager.create_a_distributed_nodal_load_by_id(loaded_nodes_v, std::set<int>{2}, std::vector<real>{input_options.get_nodal_load()}, model.glob_mesh);

    // Records
    model.scribe.track_distributed_nodes_by_id(rank, loaded_nodes_v, std::set<int>{2}, model.glob_mesh);
    time_keeper.stop_timer("bc_load_records");
    
    // Load and BC initilaisation
    time_keeper.start_timer("initialisation");
    model.initialise_restraints_n_loads();

    // initialise solution parameters 
    model.initialise_solution_parameters(input_options.max_LF, input_options.nsteps, input_options.tolerance, input_options.max_iterations);
    time_keeper.stop_timer("initialisation");
    
    // Solution
    time_keeper.start_timer("solution");
    model.solve(1);
    time_keeper.stop_timer("solution");
    time_keeper.stop_timer("all");
    // timers outputs
    #ifdef WITH_MPI
    time_keeper.log_parallel_timers(timers_names);
    model.log_parallel_timers({"U_to_nodes_mapping", 
                    "element_state_update",
                    "assembly",
                    "convergence_check",
                    "dU_calculation",
                    "material_state_update",
                    "result_recording",
                    "all"});
    #endif 
    }
}
