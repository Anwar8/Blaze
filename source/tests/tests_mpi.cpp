#include "TestHelpers.hpp"
#include "mpi.h"


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
        ASSERT_EQ(mesh.get_rank_num_nodes(), 10);
        ASSERT_EQ(mesh.get_rank_num_elems(), 9);
    }
    else if (num_ranks == 2)
    {
        ASSERT_EQ(mesh.get_rank_num_nodes(), 6);
        ASSERT_EQ(mesh.get_rank_num_elems(), 5);
    }
    else if (num_ranks == 3)
    {
        if (rank == 0)
        {
            ASSERT_EQ(mesh.get_rank_num_nodes(), 4);
            ASSERT_EQ(mesh.get_rank_num_elems(), 3);
        } else {
            ASSERT_EQ(mesh.get_rank_num_nodes(), 5);
            ASSERT_EQ(mesh.get_rank_num_elems(), 4);  
        }
    }    
    else if (num_ranks == 4)
    {
        if (rank == 0)
        {
            ASSERT_EQ(mesh.get_rank_num_nodes(), 3);
            ASSERT_EQ(mesh.get_rank_num_elems(), 2);
        } else if (rank != 3) {
            ASSERT_EQ(mesh.get_rank_num_nodes(), 4);
            ASSERT_EQ(mesh.get_rank_num_elems(), 3);  
        } else {
            ASSERT_EQ(mesh.get_rank_num_nodes(), 5);
            ASSERT_EQ(mesh.get_rank_num_elems(), 4);  
        }
    } else 
    {
        std::cout << "DistributedLineMeshTests::line_mesh_rank_counts can only run on num_ranks from 1 to 4 ranks. Got " << num_ranks << "." << std::endl;
        ASSERT_TRUE(false);
    }
}

int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  ::testing::AddGlobalTestEnvironment(new MPIEnvironment);
  return RUN_ALL_TESTS();
}