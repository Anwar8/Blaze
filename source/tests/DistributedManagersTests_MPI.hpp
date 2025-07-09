#ifndef DISTRIBUTED_MANAGERS_TESTS_HPP
#define DISTRIBUTED_MANAGERS_TESTS_HPP

#include "TestHelpers.hpp"
#define ELEMENT_TYPE NonlinearElastic


class DistributedModelFrameManagersTests : public ::testing::Test {
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
      
        }
      void TearDown() override {
  }
  };

TEST_F(DistributedModelFrameManagersTests, frame_mesh_bc_handling_count)
{   

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
    else if (num_ranks == 3)
    {
        if (rank == 0)
        {
            ASSERT_EQ(model.restraints[0].get_num_restrained_nodes(), 2);
            ASSERT_EQ(model.restraints[1].get_num_restrained_nodes(), 11);
        } 
        else if (rank == 1)
        {
            ASSERT_EQ(model.restraints[0].get_num_restrained_nodes(), 2);
            ASSERT_EQ(model.restraints[1].get_num_restrained_nodes(), 13);
            
        } 
        else if (rank == 2)
        {
            ASSERT_EQ(model.restraints[0].get_num_restrained_nodes(), 1);
            ASSERT_EQ(model.restraints[1].get_num_restrained_nodes(), 14);
        }
    }    
    else if (num_ranks == 4)
    {
        if (rank == 0)
        {
            ASSERT_EQ(model.restraints[0].get_num_restrained_nodes(), 1);
            ASSERT_EQ(model.restraints[1].get_num_restrained_nodes(), 9);
        } 
        else if (rank == 1)
        {
            ASSERT_EQ(model.restraints[0].get_num_restrained_nodes(), 1);
            ASSERT_EQ(model.restraints[1].get_num_restrained_nodes(), 11);
        } 
        else if (rank == 2)
        {
            ASSERT_EQ(model.restraints[0].get_num_restrained_nodes(), 1);
            ASSERT_EQ(model.restraints[1].get_num_restrained_nodes(), 11); 
        }
        else if (rank == 3)
        {
            ASSERT_EQ(model.restraints[0].get_num_restrained_nodes(), 1);
            ASSERT_EQ(model.restraints[1].get_num_restrained_nodes(), 9); 
        }
    }
    else if (num_ranks == 5)
    {
        if (rank == 0)
        {
            ASSERT_EQ(model.restraints[0].get_num_restrained_nodes(), 1);
            ASSERT_EQ(model.restraints[1].get_num_restrained_nodes(), 7); 
        } 
        else if (rank == 1)
        {
            ASSERT_EQ(model.restraints[0].get_num_restrained_nodes(), 1);
            ASSERT_EQ(model.restraints[1].get_num_restrained_nodes(), 10);          
        }
        else if (rank == 2)
        {
            ASSERT_EQ(model.restraints[0].get_num_restrained_nodes(), 0);
            ASSERT_EQ(model.restraints[1].get_num_restrained_nodes(), 10);          
        }
        else if (rank == 3)
        {
            ASSERT_EQ(model.restraints[0].get_num_restrained_nodes(), 1);
            ASSERT_EQ(model.restraints[1].get_num_restrained_nodes(), 9);           
        }
        else if (rank == 4)
        {
            ASSERT_EQ(model.restraints[0].get_num_restrained_nodes(), 1);
            ASSERT_EQ(model.restraints[1].get_num_restrained_nodes(), 9); 
        } 
    }
    else 
    {
        std::cout << "DistributedModelFrameManagersTests::frame_mesh_bc_handling_count can only run on num_ranks from 1 to 5 ranks. Got " << num_ranks << "." << std::endl;
        ASSERT_TRUE(false);
    }
}

  
#endif 