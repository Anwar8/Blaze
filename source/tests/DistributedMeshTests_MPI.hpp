#ifndef DISTRIBUTED_MESH_TESTS_MPI
#define DISTRIBUTED_MESH_TESTS_MPI
#include "TestHelpers.hpp"
#include "mpi.h"
/**
 * @brief got this test environment from https://bbanerjee.github.io/ParSim/mpi/c++/mpi-unit-testing-googletests-cmake/
 * @details This environment is a direct copy from Banerjee (2017). Full reference below:
 * 
 * Banerjee, Biswajit. (2017). "Unit testing with MPI, googletest, and cmake". Accessed on: 20 January 2025. URL: https://bbanerjee.github.io/ParSim/mpi/c++/mpi-unit-testing-googletests-cmake/
 */
class MPIEnvironment : public ::testing::Environment
{
public:
  virtual void SetUp() {
    char** argv;
    int argc = 0;
    int mpiError = MPI_Init(&argc, &argv);
    ASSERT_FALSE(mpiError);
  }
  virtual void TearDown() {
    int mpiError = MPI_Finalize();
    ASSERT_FALSE(mpiError);
  }
  virtual ~MPIEnvironment() {}
};

class DistributedLineMeshTests : public ::testing::Test {
  public:
    GlobalMesh mesh;
    int divisions = 9;
    std::vector<coords> end_coords = {{0.0, 0.0, 0.0}, {9.0, 0.0, 0.0}};
    std::pair<NodeIdCoordsPairsVector, ElemIdNodeIdPairVector> mesh_maps;

    void SetUp() override {
        BasicSection sect(2.06e11, 0.0125, 0.0004570000);

        mesh.set_elem_type(NonlinearElastic);
        mesh.set_basic_section(sect);
        
        mesh_maps = mesh.map_a_line_mesh(divisions, end_coords);
    }
    void TearDown() override {
}
};

TEST_F(DistributedLineMeshTests, line_mesh_rank_counts)
{   
    int rank, num_ranks;
    get_my_rank(rank);
    get_num_ranks(num_ranks);
    mesh.setup_distributed_mesh(mesh_maps.first, mesh_maps.second, rank, num_ranks);
    ASSERT_EQ(mesh.get_num_nodes(), 10);
    if (num_ranks == 1)
    {
        ASSERT_EQ(mesh.count_nodes_vector(), 10);
        ASSERT_EQ(mesh.count_elem_vector(), 9);
        ASSERT_EQ(mesh.count_interface_nodes_vector(), 0);
    }
    else if (num_ranks == 2)
    {
        ASSERT_EQ(mesh.count_nodes_vector(), 5);
        ASSERT_EQ(mesh.count_elem_vector(), 5);
        ASSERT_EQ(mesh.count_interface_nodes_vector(), 1);
    }
    else if (num_ranks == 3)
    {
        if (rank == 0)
        {
            ASSERT_EQ(mesh.count_nodes_vector(), 3);
            ASSERT_EQ(mesh.count_elem_vector(), 3);
            ASSERT_EQ(mesh.count_interface_nodes_vector(), 1);
        } 
        else if (rank == 1)
        {
            ASSERT_EQ(mesh.count_nodes_vector(), 3);
            ASSERT_EQ(mesh.count_elem_vector(), 4);
            ASSERT_EQ(mesh.count_interface_nodes_vector(), 2);
        } 
        else if (rank == 2)
        {
            ASSERT_EQ(mesh.count_nodes_vector(), 4);
            ASSERT_EQ(mesh.count_elem_vector(), 4);
            ASSERT_EQ(mesh.count_interface_nodes_vector(), 1);  
        }
    }    
    else if (num_ranks == 4)
    {
        if (rank == 0)
        {
            ASSERT_EQ(mesh.count_nodes_vector(), 2);
            ASSERT_EQ(mesh.count_elem_vector(), 2);
            ASSERT_EQ(mesh.count_interface_nodes_vector(), 1);
        } 
        else if (rank == 1 || rank == 2)
        {
            ASSERT_EQ(mesh.count_nodes_vector(), 2);
            ASSERT_EQ(mesh.count_elem_vector(), 3);
            ASSERT_EQ(mesh.count_interface_nodes_vector(), 2);
        } 
        else if (rank == 3)
        {
            ASSERT_EQ(mesh.count_nodes_vector(), 4);
            ASSERT_EQ(mesh.count_elem_vector(), 4);
            ASSERT_EQ(mesh.count_interface_nodes_vector(), 1);  
        }
    }
    else if (num_ranks == 5)
    {
        if (rank == 0 || rank == 4)
        {
            ASSERT_EQ(mesh.count_nodes_vector(), 2);
            ASSERT_EQ(mesh.count_elem_vector(), 2);
            ASSERT_EQ(mesh.count_interface_nodes_vector(), 1);
        } 
        else if (rank == 1 || rank == 2 || rank == 3)
        {
            ASSERT_EQ(mesh.count_nodes_vector(), 2);
            ASSERT_EQ(mesh.count_elem_vector(), 3);
            ASSERT_EQ(mesh.count_interface_nodes_vector(), 2);
        } 
    }
    else 
    {
        std::cout << "DistributedLineMeshTests::line_mesh_rank_counts can only run on num_ranks from 1 to 5 ranks. Got " << num_ranks << "." << std::endl;
        ASSERT_TRUE(false);
    }
}

TEST_F(DistributedLineMeshTests, line_mesh_rank_ndof_counts)
{   
    int rank, num_ranks;
    get_my_rank(rank);
    get_num_ranks(num_ranks);
    mesh.setup_distributed_mesh(mesh_maps.first, mesh_maps.second, rank, num_ranks);
    std::vector<int> ranks_ndofs_vec = mesh.get_ranks_ndofs_vector();

    if (num_ranks == 1)
    {
        check_vector_contents<int, std::vector<int>>(ranks_ndofs_vec, {60});
    }
    else if (num_ranks == 2)
    {

        check_vector_contents<int, std::vector<int>>(ranks_ndofs_vec, {30, 30});
    }
    else if (num_ranks == 3)
    {
        check_vector_contents<int, std::vector<int>>(ranks_ndofs_vec, {18, 18, 24});
    }
    else if (num_ranks == 4)
    {
        check_vector_contents<int, std::vector<int>>(ranks_ndofs_vec, {12, 12, 12, 24});
    }
    else if (num_ranks == 5)
    {
        check_vector_contents<int, std::vector<int>>(ranks_ndofs_vec, {12, 12, 12, 12, 12});
    }
    else 
    {
        std::cout << "DistributedLineMeshTests::line_mesh_rank_ndof_counts can only run on num_ranks from 1 to 5 ranks. Got " << num_ranks << "." << std::endl;
        ASSERT_TRUE(false);
    }
}

TEST_F(DistributedLineMeshTests, line_mesh_rank_interface_nzi)
{   
    int rank, num_ranks;
    get_my_rank(rank);
    get_num_ranks(num_ranks);
    mesh.setup_distributed_mesh(mesh_maps.first, mesh_maps.second, rank, num_ranks);

    if (num_ranks == 1)
    {
        ASSERT_EQ(mesh.get_node_by_record_id(1, "all")->get_nz_i(),0);
    }
    else if (num_ranks == 2)
    {
        if (rank == 0)
        {
            ASSERT_EQ(mesh.get_node_by_record_id(6, "interface")->get_nz_i(),30);
        }
        else if (rank == 1)
        {
            ASSERT_EQ(mesh.get_node_by_record_id(5, "interface")->get_nz_i(),24);
        }
    }
    else if (num_ranks == 3)
    {
        if (rank == 0)
        {
            ASSERT_EQ(mesh.get_node_by_record_id(4, "interface")->get_nz_i(),18);
        } 
        else if (rank == 1)
        {
            ASSERT_EQ(mesh.get_node_by_record_id(3, "interface")->get_nz_i(),12);
            ASSERT_EQ(mesh.get_node_by_record_id(7, "interface")->get_nz_i(),36);            
        } 
        else if (rank == 2)
        {
            ASSERT_EQ(mesh.get_node_by_record_id(6, "interface")->get_nz_i(),30);            
        }
    }    
    else if (num_ranks == 4)
    {
        if (rank == 0)
        {
            ASSERT_EQ(mesh.get_node_by_record_id(3, "interface")->get_nz_i(),12);
        } 
        else if (rank == 1)
        {
            ASSERT_EQ(mesh.get_node_by_record_id(2, "interface")->get_nz_i(),6);
            ASSERT_EQ(mesh.get_node_by_record_id(5, "interface")->get_nz_i(),24);            
        } 
        else if (rank == 2)
        {
            ASSERT_EQ(mesh.get_node_by_record_id(4, "interface")->get_nz_i(),18);
            ASSERT_EQ(mesh.get_node_by_record_id(7, "interface")->get_nz_i(),36);            
        }
        else if (rank == 3)
        {
            ASSERT_EQ(mesh.get_node_by_record_id(6, "interface")->get_nz_i(),30);           
        }
    }
    else if (num_ranks == 5)
    {
        if (rank == 0)
        {
            ASSERT_EQ(mesh.get_node_by_record_id(3, "interface")->get_nz_i(),12);
        } 
        else if (rank == 1)
        {
            ASSERT_EQ(mesh.get_node_by_record_id(2, "interface")->get_nz_i(),6);
            ASSERT_EQ(mesh.get_node_by_record_id(5, "interface")->get_nz_i(),24);            
        } 
        else if (rank == 2)
        {
            ASSERT_EQ(mesh.get_node_by_record_id(4, "interface")->get_nz_i(),18);
            ASSERT_EQ(mesh.get_node_by_record_id(7, "interface")->get_nz_i(),36);            
        }
        else if (rank == 3)
        {
            ASSERT_EQ(mesh.get_node_by_record_id(6, "interface")->get_nz_i(),30);
            ASSERT_EQ(mesh.get_node_by_record_id(9, "interface")->get_nz_i(),48);            
        }
        else if (rank == 4)
        {
            ASSERT_EQ(mesh.get_node_by_record_id(8, "interface")->get_nz_i(),42);           
        }
    }
    else 
    {
        std::cout << "DistributedLineMeshTests::line_mesh_rank_interface_nzi can only run on num_ranks from 1 to 5 ranks. Got " << num_ranks << "." << std::endl;
        ASSERT_TRUE(false);
    }
}
// TEST_F(DistributedLineMeshTests, line_mesh_rank_contents)
// {   
//     int rank, num_ranks;
//     get_my_rank(rank);
//     get_num_ranks(num_ranks);
//     mesh.setup_distributed_mesh(mesh_maps.first, mesh_maps.second, rank, num_ranks);
//     ASSERT_EQ(mesh.get_num_nodes(), 10);
//     if (num_ranks == 1)
//     {
//         ASSERT_TRUE(mesh.contains_nodes(std::set<size_t>{1, 2, 3, 4, 5, 6, 7, 8, 9, 10}));
//         ASSERT_TRUE(mesh.contains_elements(std::set<size_t>{1, 2, 3, 4, 5, 6, 7, 8, 9}));
//     }
//     else if (num_ranks == 2)
//     {
//         if (rank == 0)
//         {
//             ASSERT_TRUE(mesh.contains_nodes(std::set<size_t>{1, 2, 3, 4, 5, 6}));
//             ASSERT_TRUE(mesh.contains_elements(std::set<size_t>{1, 2, 3, 4, 5}));   
//         } else {
//             ASSERT_TRUE(mesh.contains_nodes(std::set<size_t>{5, 6, 7, 8, 9, 10}));
//             ASSERT_TRUE(mesh.contains_elements(std::set<size_t>{5, 6, 7, 8, 9}));
//         }
//     }
//     else if (num_ranks == 3)
//     {
//         if (rank == 0)
//         {
//             ASSERT_TRUE(mesh.contains_nodes(std::set<size_t>{1, 2, 3, 4}));
//             ASSERT_TRUE(mesh.contains_elements(std::set<size_t>{1, 2, 3}));
 
//         } else if (rank == 1) {
//             ASSERT_TRUE(mesh.contains_nodes(std::set<size_t>{3, 4, 5, 6, 7}));
//             ASSERT_TRUE(mesh.contains_elements(std::set<size_t>{3, 4, 5, 6}));
//         } else { 
//             ASSERT_TRUE(mesh.contains_nodes(std::set<size_t>{6, 7, 8, 9, 10}));
//             ASSERT_TRUE(mesh.contains_elements(std::set<size_t>{6, 7, 8, 9}));
//         }
//     }    
//     else if (num_ranks == 4)
//     {
//         if (rank == 0)
//         {
//             ASSERT_TRUE(mesh.contains_nodes(std::set<size_t>{1, 2, 3}));
//             ASSERT_TRUE(mesh.contains_elements(std::set<size_t>{1, 2}));
//         } else if (rank == 1) {
//             ASSERT_TRUE(mesh.contains_nodes(std::set<size_t>{2, 3, 4, 5}));
//             ASSERT_TRUE(mesh.contains_elements(std::set<size_t>{2, 3, 4}));
//         } else if (rank == 2) {
//             ASSERT_TRUE(mesh.contains_nodes(std::set<size_t>{4, 5, 6, 7}));
//             ASSERT_TRUE(mesh.contains_elements(std::set<size_t>{4, 5, 6}));
//         } else  if (rank == 3) {
//             ASSERT_TRUE(mesh.contains_nodes(std::set<size_t>{6, 7, 8, 9, 10}));
//             ASSERT_TRUE(mesh.contains_elements(std::set<size_t>{6, 7, 8, 9}));
//         }
//     } else 
//     {
//         std::cout << "DistributedLineMeshTests::line_mesh_rank_counts can only run on num_ranks from 1 to 4 ranks. Got " << num_ranks << "." << std::endl;
//         ASSERT_TRUE(false);
//     }
// }

// TEST_F(DistributedLineMeshTests, line_mesh_rank_ndofs)
// {   
//     int rank, num_ranks;
//     get_my_rank(rank);
//     get_num_ranks(num_ranks);
//     mesh.setup_distributed_mesh(mesh_maps.first, mesh_maps.second, rank, num_ranks);
//     ASSERT_EQ(mesh.get_num_nodes(), 10);
//     if (num_ranks == 1)
//     {
//         ASSERT_EQ(mesh.get_rank_ndofs(), 60);
//         ASSERT_EQ(mesh.get_rank_starting_nz_i(), 0);
//     }
//     else if (num_ranks == 2)
//     {
//         if (rank == 0)
//         {
//             ASSERT_EQ(mesh.get_rank_ndofs(), 36);
//             ASSERT_EQ(mesh.get_rank_starting_nz_i(), 0);            
//         } else {
//             ASSERT_EQ(mesh.get_rank_ndofs(), 36);
//             ASSERT_EQ(mesh.get_rank_starting_nz_i(), 30);
//         }
//     }
//     else if (num_ranks == 3)
//     {
//         if (rank == 0)
//         {
//             ASSERT_EQ(mesh.count_nodes_vector(), 4);
//             ASSERT_EQ(mesh.count_elem_vector(), 3);
//         } else {
//             ASSERT_EQ(mesh.count_nodes_vector(), 5);
//             ASSERT_EQ(mesh.count_elem_vector(), 4);  
//         }
//     }    
//     else if (num_ranks == 4)
//     {
//         if (rank == 0)
//         {
//             ASSERT_EQ(mesh.count_nodes_vector(), 3);
//             ASSERT_EQ(mesh.count_elem_vector(), 2);
//         } else if (rank != 3) {
//             ASSERT_EQ(mesh.count_nodes_vector(), 4);
//             ASSERT_EQ(mesh.count_elem_vector(), 3);  
//         } else {
//             ASSERT_EQ(mesh.count_nodes_vector(), 5);
//             ASSERT_EQ(mesh.count_elem_vector(), 4);  
//         }
//     } else 
//     {
//         std::cout << "DistributedLineMeshTests::line_mesh_rank_counts can only run on num_ranks from 1 to 4 ranks. Got " << num_ranks << "." << std::endl;
//         ASSERT_TRUE(false);
//     }
// }

#endif
