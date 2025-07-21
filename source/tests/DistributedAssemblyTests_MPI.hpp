#ifndef DISTRIBUTED_ASSEMBLY_TESTS_HPP
#define DISTRIBUTED_ASSEMBLY_TESTS_HPP

#include "TestHelpers.hpp"
#define ELEMENT_TYPE NonlinearElastic


/**
 * @brief 
 * @details this model has 32 nodes, out of which 4 have 0 DoFs, while the remaining 28 each have 3 DoFs for a total of 84 DoFs in the whole system. 
 * 
 */
class DistributedModelFrameAssemblyTests : public ::testing::Test {
    public:
      Model model;

      int rank = 0;
      int num_ranks = 0;
      int nbays = 3;
      int nfloors = 2;
      real bay_length = 6;
      real floor_height = 4;
      int beam_divisions = 3;
      int column_divisions = 2;
      CommonSectionDefinitions common;

      FrameMesh the_frame;
      NodalRestraint column_bases;
      NodalRestraint out_of_plane_restraint; 
      
      void SetUp() override {
          common.initialise_section();
          get_my_rank(rank);
          get_num_ranks(num_ranks);
          model.create_distributed_frame_mesh(nbays, nfloors, bay_length, floor_height, beam_divisions, column_divisions, NonlinearPlastic, common.I_section);
          

          the_frame = model.glob_mesh.get_frame();
          

            column_bases.assign_dofs_restraints(std::set<int>{0, 1, 2, 3, 4, 5}); // fixed support
            column_bases.assign_distributed_nodes_by_record_id(the_frame.get_column_bases(), model.glob_mesh);

            NodalRestraint out_of_plane_restraint;  
            out_of_plane_restraint.assign_dofs_restraints(std::set<int>{1, 3, 4});
            out_of_plane_restraint.assign_distributed_nodes_by_record_id(the_frame.get_out_of_plane_nodes(), model.glob_mesh);

            model.restraints.push_back(column_bases);
            model.restraints.push_back(out_of_plane_restraint);

            std::set<unsigned> loaded_nodes = the_frame.get_all_beam_line_node_ids(true);
            std::vector<unsigned> loaded_nodes_v = std::vector<unsigned>(loaded_nodes.begin(), loaded_nodes.end());
            model.load_manager.create_a_distributed_nodal_load_by_id(loaded_nodes_v, std::set<int>{2}, std::vector<real>{-1000}, model.glob_mesh);
            model.initialise_restraints_n_loads();
            
        }
      void TearDown() override {
  }
  };

TEST_F(DistributedModelFrameAssemblyTests, frame_mesh_vector_sizes)
{   
    // std::cout << std::endl << std::endl;
    // std::cout << "The fully-assembled stiffness matrix K is:" << std::endl;
    // model.assembler.print_distributed_maths_object("K", Teuchos::VERB_EXTREME);

    // // initialise solution parameters
    real max_LF = 1;
    int nsteps = 100;
    real tolerance = 1e-2;
    int max_iterations = 10;
    model.initialise_solution_parameters(max_LF, nsteps, tolerance, max_iterations);
    model.solve(1);
    std::cout << "SECTION:Timing" << std::endl;
    model.log_parallel_timers({"U_to_nodes_mapping", 
                    "element_state_update", 
                    "element_global_response",
                    "assembly",
                    "convergence_check",
                    "dU_calculation",
                    "material_state_update",
                    "result_recording",
                    "all"});


    if (num_ranks == 1)
    {
        ASSERT_EQ(model.restraints[0].get_num_restrained_nodes(), 4);
        ASSERT_EQ(model.restraints[1].get_num_restrained_nodes(), 28);
    }
    else if (num_ranks == 2)
    {

        ASSERT_EQ(model.restraints[0].get_num_restrained_nodes(), 2);
        // remember that the interface nodes are also restrained, so we have two more nodes 
        ASSERT_EQ(model.restraints[1].get_num_restrained_nodes(), 16);

    }
    else 
    {
        std::cout << "DistributedModelFrameAssemblyTests::frame_mesh_vector_sizes can only run on num_ranks from 1 to 5 ranks. Got " << num_ranks << "." << std::endl;
        ASSERT_TRUE(false);
    }
}
  

class DistributedModelSimplySuportedUdlPlastic : public ::testing::Test {
    public:
      Model model;
      CommonSectionDefinitions common;
      int divisions = 100;
      real beam_length = 5.0;
      real y_udl = -1e4; // N/m
      
      unsigned mid_node = (divisions/2) + 1;
      std::vector<unsigned> loaded_nodes = std::vector<unsigned>(divisions - 1);
      std::set<unsigned> end_nodes = {1, unsigned(1 + divisions)};
      int tracked_dof = 2;
      
      int rank = 0;
      int num_ranks = 0;
      
      NodalRestraint end_restraints_1;
      NodalRestraint end_restraints_2;
      NodalRestraint out_of_plane_restraint; 
      
      void SetUp() override {
          common.initialise_section();
          get_my_rank(rank);
          get_num_ranks(num_ranks);
          model.create_distributed_line_mesh(divisions, {{0.0, 0.0, 0.0}, {beam_length, 0.0, 0.0}}, NonlinearPlastic, common.I_section);

          std::set<unsigned> end_nodes_1 = {1};
          std::set<unsigned> end_nodes_2 = {unsigned(1 + divisions)};
          end_restraints_1.assign_dofs_restraints(std::set<int>{0, 1, 2, 3, 4}); // restrain x translation, x rotation, y translation, y rotation, and z translation
          end_restraints_1.assign_distributed_nodes_by_record_id(end_nodes_1, model.glob_mesh);
          
          end_restraints_2.assign_dofs_restraints(std::set<int>{1, 2, 3, 4}); // restrain x rotation, y translation, y 
          end_restraints_2.assign_distributed_nodes_by_record_id(end_nodes_2, model.glob_mesh);

          // create the loaded and restrained intermediate nodes
          std::iota(loaded_nodes.begin(), loaded_nodes.end(), 2);
          out_of_plane_restraint.assign_dofs_restraints(std::set<int>{1, 3, 4}); // restrain x rotation, y rotation, and z translation
          out_of_plane_restraint.assign_distributed_nodes_by_record_id(loaded_nodes, model.glob_mesh);

          model.restraints.push_back(end_restraints_1);
          model.restraints.push_back(end_restraints_2);
          model.restraints.push_back(out_of_plane_restraint);

        // calculate and apply load
          real y_load = y_udl*beam_length/(divisions - 1);
          model.load_manager.create_a_distributed_nodal_load_by_id(loaded_nodes, std::set<int>{tracked_dof}, std::vector<real>{y_load}, model.glob_mesh);
        
        // track mid-point
          model.scribe.track_distributed_nodes_by_id(rank ,std::vector<unsigned>{mid_node}, std::set<int>{tracked_dof}, model.glob_mesh);
        // initialise solver and solve
          model.initialise_restraints_n_loads();
          model.initialise_solution_parameters(1.0, 100, 1e-3, 10);
          model.solve(-1);
            
        }
      void TearDown() override {
  }
  };

  TEST_F(DistributedModelSimplySuportedUdlPlastic, CheckMidSpanDeflection)
  {
    if (model.glob_mesh.owns_node_record_id(mid_node))
    {
        std::shared_ptr<Node> node = model.glob_mesh.get_node_by_id(mid_node);
  
        std::vector<Record> record_library = model.scribe.get_record_library();
        Record record = record_library.back();
    
        std::array<std::vector<real>, 6> recorded_data = record.get_recorded_data(); 
        std::vector<real> disp_data = recorded_data[tracked_dof];

        // $\delta = \frac{5 w L^4}{384 EI} $
        real correct_disp = 5*y_udl*std::pow(beam_length, 4)/(384*(YOUNGS_MODULUS)*(common.moment_of_inertia));
        real tolerance = std::abs(PERCENT_TOLERANCE*correct_disp);
        EXPECT_NEAR(disp_data.back(), correct_disp, tolerance);
    } else {
        GTEST_SKIP() << "Rank " << rank  << " does not own the node; skipping check.";
    }

  }


  class DistributedModelSimplySuportedUdlElastic : public ::testing::Test {
    public:
      Model model;
      CommonSectionDefinitions common;
      
      int divisions = 100;
      real beam_length = 5.0;
      real y_udl = -1e4; // N/m
      
      unsigned mid_node = (divisions/2) + 1;
      std::vector<unsigned> loaded_nodes = std::vector<unsigned>(divisions - 1);
      std::set<unsigned> end_nodes = {1, unsigned(1 + divisions)};
      int tracked_dof = 2;
      
      int rank = 0;
      int num_ranks = 0;
      
      NodalRestraint end_restraints_1;
      NodalRestraint end_restraints_2;
      NodalRestraint out_of_plane_restraint; 
      
      void SetUp() override {
          common.initialise_section();
          get_my_rank(rank);
          get_num_ranks(num_ranks);
        //   model.create_distributed_line_mesh(divisions, {{0.0, 0.0, 0.0}, {beam_length, 0.0, 0.0}}, NonlinearPlastic, common.I_section);
          BasicSection sect(2.06e11, 0.0125, 0.0004570000);
          model.create_distributed_line_mesh(divisions, {{0.0, 0.0, 0.0}, {beam_length, 0.0, 0.0}}, LinearElastic, sect);

          std::set<unsigned> end_nodes_1 = {1};
          std::set<unsigned> end_nodes_2 = {unsigned(1 + divisions)};
          end_restraints_1.assign_dofs_restraints(std::set<int>{0, 1, 2, 3, 4}); // restrain x translation, x rotation, y translation, y rotation, and z translation
          end_restraints_1.assign_distributed_nodes_by_record_id(end_nodes_1, model.glob_mesh);
          
          end_restraints_2.assign_dofs_restraints(std::set<int>{1, 2, 3, 4}); // restrain x rotation, y translation, y 
          end_restraints_2.assign_distributed_nodes_by_record_id(end_nodes_2, model.glob_mesh);

          // create the loaded and restrained intermediate nodes
          std::iota(loaded_nodes.begin(), loaded_nodes.end(), 2);
          out_of_plane_restraint.assign_dofs_restraints(std::set<int>{1, 3, 4}); // restrain x rotation, y rotation, and z translation
          out_of_plane_restraint.assign_distributed_nodes_by_record_id(loaded_nodes, model.glob_mesh);

          model.restraints.push_back(end_restraints_1);
          model.restraints.push_back(end_restraints_2);
          model.restraints.push_back(out_of_plane_restraint);

        // calculate and apply load
          real y_load = y_udl*beam_length/(divisions - 1);
          model.load_manager.create_a_distributed_nodal_load_by_id(loaded_nodes, std::set<int>{tracked_dof}, std::vector<real>{y_load}, model.glob_mesh);
        
        // track mid-point
          model.scribe.track_distributed_nodes_by_id(rank ,std::vector<unsigned>{mid_node}, std::set<int>{tracked_dof}, model.glob_mesh);
        // initialise solver and solve
          model.initialise_restraints_n_loads();
          model.initialise_solution_parameters(1.0, 100, 1e-3, 10);
          model.solve(-1);
            
        }
      void TearDown() override {
  }
  };

  TEST_F(DistributedModelSimplySuportedUdlElastic, CheckMidSpanDeflection)
  {
    if (model.glob_mesh.owns_node_record_id(mid_node))
    {
        std::shared_ptr<Node> node = model.glob_mesh.get_node_by_id(mid_node);
  
        std::vector<Record> record_library = model.scribe.get_record_library();
        Record record = record_library.back();
    
        std::array<std::vector<real>, 6> recorded_data = record.get_recorded_data(); 
        std::vector<real> disp_data = recorded_data[tracked_dof];

        // $\delta = \frac{5 w L^4}{384 EI} $
        real correct_disp = 5*y_udl*std::pow(beam_length, 4)/(384*(YOUNGS_MODULUS)*(common.moment_of_inertia));
        real tolerance = std::abs(PERCENT_TOLERANCE*correct_disp);
        EXPECT_NEAR(disp_data.back(), correct_disp, tolerance);
    } else {
        GTEST_SKIP() << "Rank " << rank  << " does not own the node; skipping check.";
    }

  }
#endif 