#ifndef DISTRIBUTED_MESH_TESTS_HPP
#define DISTRIBUTED_MESH_TESTS_HPP

#include "TestHelpers.hpp"
#define ELEMENT_TYPE NonlinearElastic


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

/**
 * Start testing all the helper functions not necessarily the actual GlobalMesh object itself, yet.
 * The MPI_Allgather call should work find if we compile with MPI and just run normally as it would have a num_ranks = 1, and rank = 0.
 */

/**
 * @name populate_node_rank_maps
 * @brief Tests the populate_node_rank_maps function with different num_ranks.
 */
//@{
TEST_F(DistributedLineMeshTests, populate_node_rank_maps_10_on_1_counts)
{
    std::map<size_t, int> node_rank_map; 
    std::map<int, std::set<size_t>> rank_nodes_map;
    mesh.populate_node_rank_maps(node_rank_map, rank_nodes_map, mesh_maps.first, 1);
    EXPECT_EQ(rank_nodes_map.size(), 1);
    EXPECT_EQ(node_rank_map.size(), 10);
}

TEST_F(DistributedLineMeshTests, populate_node_rank_maps_10_on_1_values)
{
    std::map<size_t, int> node_rank_map; 
    std::map<int, std::set<size_t>> rank_nodes_map;
    mesh.populate_node_rank_maps(node_rank_map, rank_nodes_map, mesh_maps.first, 1);

    check_set_ids(rank_nodes_map[0], 0, 10, 1);
    
    for (int i = 1; i <= 10; ++i)
    {
        EXPECT_EQ(node_rank_map[i], 0);
    }
}

TEST_F(DistributedLineMeshTests, populate_node_rank_maps_10_on_2_counts)
{
    std::map<size_t, int> node_rank_map; 
    std::map<int, std::set<size_t>> rank_nodes_map;
    mesh.populate_node_rank_maps(node_rank_map, rank_nodes_map, mesh_maps.first, 2);
    EXPECT_EQ(rank_nodes_map.size(), 2);
    EXPECT_EQ(node_rank_map.size(), 10);
}

TEST_F(DistributedLineMeshTests, populate_node_rank_maps_10_on_2_values)
{
    std::map<size_t, int> node_rank_map; 
    std::map<int, std::set<size_t>> rank_nodes_map;
    mesh.populate_node_rank_maps(node_rank_map, rank_nodes_map, mesh_maps.first, 2);

    check_set_ids(rank_nodes_map[0], 0, 5, 1);
    check_set_ids(rank_nodes_map[1], 0, 5, 6);
    
    for (int i = 1; i <= 10; ++i)
    {
        if (i <= 5)
        {
            EXPECT_EQ(node_rank_map[i], 0);
        } else {
            EXPECT_EQ(node_rank_map[i], 1);
        }
    }
}

TEST_F(DistributedLineMeshTests, populate_node_rank_maps_10_on_3_counts)
{
    std::map<size_t, int> node_rank_map; 
    std::map<int, std::set<size_t>> rank_nodes_map;
    mesh.populate_node_rank_maps(node_rank_map, rank_nodes_map, mesh_maps.first, 3);
    EXPECT_EQ(rank_nodes_map.size(), 3);
    EXPECT_EQ(node_rank_map.size(), 10);
}

TEST_F(DistributedLineMeshTests, populate_node_rank_maps_10_on_3_values)
{
    std::map<size_t, int> node_rank_map; 
    std::map<int, std::set<size_t>> rank_nodes_map;
    mesh.populate_node_rank_maps(node_rank_map, rank_nodes_map, mesh_maps.first, 3);

    check_set_ids(rank_nodes_map[0], 0, 3, 1);
    check_set_ids(rank_nodes_map[1], 0, 3, 4);
    check_set_ids(rank_nodes_map[2], 0, 4, 7);
    
    for (int i = 1; i <= 10; ++i)
    {
        if (i <= 3)
        {
            EXPECT_EQ(node_rank_map[i], 0);
        } 
        else if (i > 3 && i <= 6) 
        {
            EXPECT_EQ(node_rank_map[i], 1);
        } 
        else 
        {
            EXPECT_EQ(node_rank_map[i], 2);
        }
    }
}
//@}
/**
 * @name populate_node_element_map
 * @brief Tests the populate_node_element_map function. No variations needed.
 */
//@{

TEST_F(DistributedLineMeshTests, populate_node_element_map_9_elems_counts)
{
    std::map<size_t, std::set<size_t>> node_element_map;
    mesh.populate_node_element_map(node_element_map, mesh_maps.second);
    EXPECT_EQ(node_element_map.size(), 10);
    EXPECT_EQ(node_element_map[1].size(), 1);
    EXPECT_EQ(node_element_map[10].size(), 1);
    for (int i = 2; i < 10; ++i)
    {
        EXPECT_EQ(node_element_map[i].size(), 2);
    }
}

TEST_F(DistributedLineMeshTests, populate_node_element_map_9_elems_values)
{
    std::map<size_t, std::set<size_t>> node_element_map;
    mesh.populate_node_element_map(node_element_map, mesh_maps.second);
    EXPECT_EQ(*node_element_map[1].begin(), 1);
    EXPECT_EQ(*node_element_map[10].begin(), 9);
    
    for (int i = 2; i < 10; ++i)
    {
        auto elem_set_it = node_element_map[i].begin();
        EXPECT_EQ(*elem_set_it, i - 1);
        elem_set_it++;
        EXPECT_EQ(*elem_set_it, i);
    }
}
//@}
/**
 * @name find_rank_elements
 * @brief Tests the find_rank_elements function with different num_ranks and for each rank. Must have the correct number of elements and element ids for each rank.
 * @details this is a unique and difficult to understand test. This because of the way distribution works in Blaze. Each rank "owns" nodes, not elements, and this ownership governs which elements will be created on which ranks which this function aims to study. Therefore, for the case when the number of nodes is indivisible by number of ranks, the number of elements on each rank will vary and be quite unintuitive. For example, 10 nodes on 3 ranks will results in nodes (1, 2, 3) and elements (1, 2, 3) on rank 0, nodes (4, 5, 6) and elements (3, 4, 5, 6) on rank 1, and nodes (7, 8, 9, 10) and elements (6, 7, 8, 9) on rank 2. It may be required to draw this to understand it.
 */
//@{
TEST_F(DistributedLineMeshTests, find_rank_elements_10_on_1_counts)
{
    std::map<size_t, int> node_rank_map; 
    std::map<int, std::set<size_t>> rank_nodes_map;
    mesh.populate_node_rank_maps(node_rank_map, rank_nodes_map, mesh_maps.first, 1);
    std::map<size_t, std::set<size_t>> node_element_map;
    mesh.populate_node_element_map(node_element_map, mesh_maps.second);

    std::set<size_t> elem_id_set_on_rank_0;
    mesh.find_rank_elements(elem_id_set_on_rank_0, rank_nodes_map, node_element_map, 0);
    int rank_nelems = elem_id_set_on_rank_0.size();

    EXPECT_EQ(rank_nelems, 9);
}

TEST_F(DistributedLineMeshTests, find_rank_elements_10_on_1_values)
{
    std::map<size_t, int> node_rank_map; 
    std::map<int, std::set<size_t>> rank_nodes_map;
    mesh.populate_node_rank_maps(node_rank_map, rank_nodes_map, mesh_maps.first, 1);
    std::map<size_t, std::set<size_t>> node_element_map;
    mesh.populate_node_element_map(node_element_map, mesh_maps.second);

    std::set<size_t> elem_id_set_on_rank_0;
    mesh.find_rank_elements(elem_id_set_on_rank_0, rank_nodes_map, node_element_map, 0);

    check_set_ids(elem_id_set_on_rank_0, 0, 9, 1);
}

TEST_F(DistributedLineMeshTests, find_rank_elements_10_on_2_counts)
{
    std::map<size_t, int> node_rank_map; 
    std::map<int, std::set<size_t>> rank_nodes_map;
    mesh.populate_node_rank_maps(node_rank_map, rank_nodes_map, mesh_maps.first, 2);
    std::map<size_t, std::set<size_t>> node_element_map;
    mesh.populate_node_element_map(node_element_map, mesh_maps.second);

    std::set<size_t> elem_id_set_on_rank_0;
    mesh.find_rank_elements(elem_id_set_on_rank_0, rank_nodes_map, node_element_map, 0);
    int rank0_nelems = elem_id_set_on_rank_0.size();

    std::set<size_t> elem_id_set_on_rank_1;
    mesh.find_rank_elements(elem_id_set_on_rank_1, rank_nodes_map, node_element_map, 1);
    int rank1_nelems = elem_id_set_on_rank_1.size();

    EXPECT_EQ(rank0_nelems, 5);
    EXPECT_EQ(rank1_nelems, 5);
}

TEST_F(DistributedLineMeshTests, find_rank_elements_10_on_2_values)
{
    std::map<size_t, int> node_rank_map; 
    std::map<int, std::set<size_t>> rank_nodes_map;
    mesh.populate_node_rank_maps(node_rank_map, rank_nodes_map, mesh_maps.first, 2);
    std::map<size_t, std::set<size_t>> node_element_map;
    mesh.populate_node_element_map(node_element_map, mesh_maps.second);

    std::set<size_t> elem_id_set_on_rank_0;
    mesh.find_rank_elements(elem_id_set_on_rank_0, rank_nodes_map, node_element_map, 0);

    std::set<size_t> elem_id_set_on_rank_1;
    mesh.find_rank_elements(elem_id_set_on_rank_1, rank_nodes_map, node_element_map, 1);

    check_set_ids(elem_id_set_on_rank_0, 0, 5, 1);
    check_set_ids(elem_id_set_on_rank_1, 0, 5, 5);
}

TEST_F(DistributedLineMeshTests, find_rank_elements_10_on_3_counts)
{
    std::map<size_t, int> node_rank_map; 
    std::map<int, std::set<size_t>> rank_nodes_map;
    mesh.populate_node_rank_maps(node_rank_map, rank_nodes_map, mesh_maps.first, 3);
    std::map<size_t, std::set<size_t>> node_element_map;
    mesh.populate_node_element_map(node_element_map, mesh_maps.second);

    std::set<size_t> elem_id_set_on_rank_0;
    mesh.find_rank_elements(elem_id_set_on_rank_0, rank_nodes_map, node_element_map, 0);
    int rank0_nelems = elem_id_set_on_rank_0.size();

    std::set<size_t> elem_id_set_on_rank_1;
    mesh.find_rank_elements(elem_id_set_on_rank_1, rank_nodes_map, node_element_map, 1);
    int rank1_nelems = elem_id_set_on_rank_1.size();

    std::set<size_t> elem_id_set_on_rank_2;
    mesh.find_rank_elements(elem_id_set_on_rank_2, rank_nodes_map, node_element_map, 2);
    int rank2_nelems = elem_id_set_on_rank_2.size();

    EXPECT_EQ(rank0_nelems, 3);
    EXPECT_EQ(rank1_nelems, 4);
    EXPECT_EQ(rank2_nelems, 4);
}

TEST_F(DistributedLineMeshTests, find_rank_elements_10_on_3_values)
{
    std::map<size_t, int> node_rank_map; 
    std::map<int, std::set<size_t>> rank_nodes_map;
    mesh.populate_node_rank_maps(node_rank_map, rank_nodes_map, mesh_maps.first, 3);
    std::map<size_t, std::set<size_t>> node_element_map;
    mesh.populate_node_element_map(node_element_map, mesh_maps.second);

    std::set<size_t> elem_id_set_on_rank_0;
    mesh.find_rank_elements(elem_id_set_on_rank_0, rank_nodes_map, node_element_map, 0);

    std::set<size_t> elem_id_set_on_rank_1;
    mesh.find_rank_elements(elem_id_set_on_rank_1, rank_nodes_map, node_element_map, 1);

    std::set<size_t> elem_id_set_on_rank_2;
    mesh.find_rank_elements(elem_id_set_on_rank_2, rank_nodes_map, node_element_map, 2);

    check_set_ids(elem_id_set_on_rank_0, 0, 3, 1);
    check_set_ids(elem_id_set_on_rank_1, 0, 4, 3);
    check_set_ids(elem_id_set_on_rank_2, 0, 4, 6);
}
//@}

/**
 * @name filter_element_vector
 * @brief Tests the filter_element_vector function with different num_ranks and for each rank. Must have the correct number of elements and element ids for each rank.
 */
//@{
TEST_F(DistributedLineMeshTests, filter_element_vector_10_on_1_counts)
{
    std::map<size_t, int> node_rank_map; 
    std::map<int, std::set<size_t>> rank_nodes_map;
    mesh.populate_node_rank_maps(node_rank_map, rank_nodes_map, mesh_maps.first, 1);
    std::map<size_t, std::set<size_t>> node_element_map;
    mesh.populate_node_element_map(node_element_map, mesh_maps.second);

    std::set<size_t> elem_id_set_on_rank_0;
    mesh.find_rank_elements(elem_id_set_on_rank_0, rank_nodes_map, node_element_map, 0);
    int rank_nelems = elem_id_set_on_rank_0.size();

    ElemIdNodeIdPairVector elem_nodes_vector_on_rank_0;
    elem_nodes_vector_on_rank_0.reserve(rank_nelems);
    mesh.filter_element_vector(elem_nodes_vector_on_rank_0, elem_id_set_on_rank_0, mesh_maps.second, 0);

    EXPECT_EQ(elem_nodes_vector_on_rank_0.size(), 9);
}

TEST_F(DistributedLineMeshTests, filter_element_vector_10_on_1_values)
{
    std::map<size_t, int> node_rank_map; 
    std::map<int, std::set<size_t>> rank_nodes_map;
    mesh.populate_node_rank_maps(node_rank_map, rank_nodes_map, mesh_maps.first, 1);
    std::map<size_t, std::set<size_t>> node_element_map;
    mesh.populate_node_element_map(node_element_map, mesh_maps.second);

    std::set<size_t> elem_id_set_on_rank_0;
    mesh.find_rank_elements(elem_id_set_on_rank_0, rank_nodes_map, node_element_map, 0);
    int rank0_nelems = elem_id_set_on_rank_0.size();

    ElemIdNodeIdPairVector elem_nodes_vector_on_rank_0;
    elem_nodes_vector_on_rank_0.reserve(rank0_nelems);
    mesh.filter_element_vector(elem_nodes_vector_on_rank_0, elem_id_set_on_rank_0, mesh_maps.second, 0);

    check_vector_ids(elem_nodes_vector_on_rank_0, 0, 9, 1);
}

TEST_F(DistributedLineMeshTests, filter_element_vector_10_on_2_counts)
{
    std::map<size_t, int> node_rank_map; 
    std::map<int, std::set<size_t>> rank_nodes_map;
    mesh.populate_node_rank_maps(node_rank_map, rank_nodes_map, mesh_maps.first, 2);
    std::map<size_t, std::set<size_t>> node_element_map;
    mesh.populate_node_element_map(node_element_map, mesh_maps.second);

    std::set<size_t> elem_id_set_on_rank_0;
    mesh.find_rank_elements(elem_id_set_on_rank_0, rank_nodes_map, node_element_map, 0);
    int rank0_nelems = elem_id_set_on_rank_0.size();

    std::set<size_t> elem_id_set_on_rank_1;
    mesh.find_rank_elements(elem_id_set_on_rank_1, rank_nodes_map, node_element_map, 1);
    int rank1_nelems = elem_id_set_on_rank_1.size();

    ElemIdNodeIdPairVector elem_nodes_vector_on_rank_0;
    elem_nodes_vector_on_rank_0.reserve(rank0_nelems);
    mesh.filter_element_vector(elem_nodes_vector_on_rank_0, elem_id_set_on_rank_0, mesh_maps.second, 0);

    ElemIdNodeIdPairVector elem_nodes_vector_on_rank_1;
    elem_nodes_vector_on_rank_1.reserve(rank1_nelems);
    mesh.filter_element_vector(elem_nodes_vector_on_rank_1, elem_id_set_on_rank_1, mesh_maps.second, 1);

    EXPECT_EQ(elem_nodes_vector_on_rank_0.size(), 5);
    EXPECT_EQ(elem_nodes_vector_on_rank_1.size(), 5);
}

TEST_F(DistributedLineMeshTests, filter_element_vector_10_on_2_values)
{
    std::map<size_t, int> node_rank_map; 
    std::map<int, std::set<size_t>> rank_nodes_map;
    mesh.populate_node_rank_maps(node_rank_map, rank_nodes_map, mesh_maps.first, 2);
    std::map<size_t, std::set<size_t>> node_element_map;
    mesh.populate_node_element_map(node_element_map, mesh_maps.second);

    std::set<size_t> elem_id_set_on_rank_0;
    mesh.find_rank_elements(elem_id_set_on_rank_0, rank_nodes_map, node_element_map, 0);
    int rank0_nelems = elem_id_set_on_rank_0.size();

    std::set<size_t> elem_id_set_on_rank_1;
    mesh.find_rank_elements(elem_id_set_on_rank_1, rank_nodes_map, node_element_map, 1);
    int rank1_nelems = elem_id_set_on_rank_1.size();

    ElemIdNodeIdPairVector elem_nodes_vector_on_rank_0;
    elem_nodes_vector_on_rank_0.reserve(rank0_nelems);
    mesh.filter_element_vector(elem_nodes_vector_on_rank_0, elem_id_set_on_rank_0, mesh_maps.second, 0);

    ElemIdNodeIdPairVector elem_nodes_vector_on_rank_1;
    elem_nodes_vector_on_rank_1.reserve(rank1_nelems);
    mesh.filter_element_vector(elem_nodes_vector_on_rank_1, elem_id_set_on_rank_1, mesh_maps.second, 1);

    check_vector_ids(elem_nodes_vector_on_rank_0, 0, 5, 1);
    check_vector_ids(elem_nodes_vector_on_rank_1, 0, 5, 5);
}

TEST_F(DistributedLineMeshTests, filter_element_vector_10_on_3_counts)
{
    std::map<size_t, int> node_rank_map; 
    std::map<int, std::set<size_t>> rank_nodes_map;
    mesh.populate_node_rank_maps(node_rank_map, rank_nodes_map, mesh_maps.first, 3);
    std::map<size_t, std::set<size_t>> node_element_map;
    mesh.populate_node_element_map(node_element_map, mesh_maps.second);

    std::set<size_t> elem_id_set_on_rank_0;
    mesh.find_rank_elements(elem_id_set_on_rank_0, rank_nodes_map, node_element_map, 0);
    int rank0_nelems = elem_id_set_on_rank_0.size();

    std::set<size_t> elem_id_set_on_rank_1;
    mesh.find_rank_elements(elem_id_set_on_rank_1, rank_nodes_map, node_element_map, 1);
    int rank1_nelems = elem_id_set_on_rank_1.size();

    std::set<size_t> elem_id_set_on_rank_2;
    mesh.find_rank_elements(elem_id_set_on_rank_2, rank_nodes_map, node_element_map, 2);
    int rank2_nelems = elem_id_set_on_rank_2.size();

    ElemIdNodeIdPairVector elem_nodes_vector_on_rank_0;
    elem_nodes_vector_on_rank_0.reserve(rank0_nelems);
    mesh.filter_element_vector(elem_nodes_vector_on_rank_0, elem_id_set_on_rank_0, mesh_maps.second, 0);

    ElemIdNodeIdPairVector elem_nodes_vector_on_rank_1;
    elem_nodes_vector_on_rank_1.reserve(rank1_nelems);
    mesh.filter_element_vector(elem_nodes_vector_on_rank_1, elem_id_set_on_rank_1, mesh_maps.second, 1);

    ElemIdNodeIdPairVector elem_nodes_vector_on_rank_2;
    elem_nodes_vector_on_rank_2.reserve(rank2_nelems);
    mesh.filter_element_vector(elem_nodes_vector_on_rank_2, elem_id_set_on_rank_2, mesh_maps.second, 2);

    EXPECT_EQ(elem_nodes_vector_on_rank_0.size(), 3);
    EXPECT_EQ(elem_nodes_vector_on_rank_1.size(), 4);
    EXPECT_EQ(elem_nodes_vector_on_rank_2.size(), 4);
}

TEST_F(DistributedLineMeshTests, filter_element_vector_10_on_3_values)
{
    std::map<size_t, int> node_rank_map; 
    std::map<int, std::set<size_t>> rank_nodes_map;
    mesh.populate_node_rank_maps(node_rank_map, rank_nodes_map, mesh_maps.first, 3);
    std::map<size_t, std::set<size_t>> node_element_map;
    mesh.populate_node_element_map(node_element_map, mesh_maps.second);

    std::set<size_t> elem_id_set_on_rank_0;
    mesh.find_rank_elements(elem_id_set_on_rank_0, rank_nodes_map, node_element_map, 0);
    int rank0_nelems = elem_id_set_on_rank_0.size();

    std::set<size_t> elem_id_set_on_rank_1;
    mesh.find_rank_elements(elem_id_set_on_rank_1, rank_nodes_map, node_element_map, 1);
    int rank1_nelems = elem_id_set_on_rank_1.size();

    std::set<size_t> elem_id_set_on_rank_2;
    mesh.find_rank_elements(elem_id_set_on_rank_2, rank_nodes_map, node_element_map, 2);
    int rank2_nelems = elem_id_set_on_rank_2.size();

    ElemIdNodeIdPairVector elem_nodes_vector_on_rank_0;
    elem_nodes_vector_on_rank_0.reserve(rank0_nelems);
    mesh.filter_element_vector(elem_nodes_vector_on_rank_0, elem_id_set_on_rank_0, mesh_maps.second, 0);

    ElemIdNodeIdPairVector elem_nodes_vector_on_rank_1;
    elem_nodes_vector_on_rank_1.reserve(rank1_nelems);
    mesh.filter_element_vector(elem_nodes_vector_on_rank_1, elem_id_set_on_rank_1, mesh_maps.second, 1);

    ElemIdNodeIdPairVector elem_nodes_vector_on_rank_2;
    elem_nodes_vector_on_rank_2.reserve(rank2_nelems);
    mesh.filter_element_vector(elem_nodes_vector_on_rank_2, elem_id_set_on_rank_2, mesh_maps.second, 2);

    check_vector_ids(elem_nodes_vector_on_rank_0, 0, 3, 1);
    check_vector_ids(elem_nodes_vector_on_rank_1, 0, 4, 3);
    check_vector_ids(elem_nodes_vector_on_rank_2, 0, 4, 6);
}
//@}

/**
 * @name find_rank_nodes
 * @brief Tests the find_rank_nodes function with different num_ranks and for each rank. Must correctly identify which nodes need to be replicated, and which nodes are retained from original. 
 * @details See the description of \ref find_rank_elements for more details on the kind of unintuitive behaviour expected as this function will need to get even more nodes to allow for creation of the replicated elements created in find_rank_elements.
 */
//@{
TEST_F(DistributedLineMeshTests, find_rank_nodes_10_on_1_counts)
{
    std::map<size_t, int> node_rank_map; 
    std::map<int, std::set<size_t>> rank_nodes_map;
    mesh.populate_node_rank_maps(node_rank_map, rank_nodes_map, mesh_maps.first, 1);
    std::map<size_t, std::set<size_t>> node_element_map;
    mesh.populate_node_element_map(node_element_map, mesh_maps.second);

    std::set<size_t> elem_id_set_on_rank_0;
    mesh.find_rank_elements(elem_id_set_on_rank_0, rank_nodes_map, node_element_map, 0);
    int rank_nelems = elem_id_set_on_rank_0.size();

    ElemIdNodeIdPairVector elem_nodes_vector_on_rank_0;
    elem_nodes_vector_on_rank_0.reserve(rank_nelems);
    mesh.filter_element_vector(elem_nodes_vector_on_rank_0, elem_id_set_on_rank_0, mesh_maps.second, 0);

    std::set<size_t> node_id_set_on_rank_0;
    mesh.find_rank_nodes(node_id_set_on_rank_0, elem_nodes_vector_on_rank_0);
    int rank0_nnodes = node_id_set_on_rank_0.size();

    EXPECT_EQ(rank0_nnodes, 10);
}

TEST_F(DistributedLineMeshTests, find_rank_nodes_10_on_1_values)
{
    std::map<size_t, int> node_rank_map; 
    std::map<int, std::set<size_t>> rank_nodes_map;
    mesh.populate_node_rank_maps(node_rank_map, rank_nodes_map, mesh_maps.first, 1);
    std::map<size_t, std::set<size_t>> node_element_map;
    mesh.populate_node_element_map(node_element_map, mesh_maps.second);

    std::set<size_t> elem_id_set_on_rank_0;
    mesh.find_rank_elements(elem_id_set_on_rank_0, rank_nodes_map, node_element_map, 0);
    int rank0_nelems = elem_id_set_on_rank_0.size();

    ElemIdNodeIdPairVector elem_nodes_vector_on_rank_0;
    elem_nodes_vector_on_rank_0.reserve(rank0_nelems);
    mesh.filter_element_vector(elem_nodes_vector_on_rank_0, elem_id_set_on_rank_0, mesh_maps.second, 0);

    std::set<size_t> node_id_set_on_rank_0;
    mesh.find_rank_nodes(node_id_set_on_rank_0, elem_nodes_vector_on_rank_0);

    check_set_ids(node_id_set_on_rank_0, 0, 10, 1);
}

TEST_F(DistributedLineMeshTests, find_rank_nodes_10_on_2_counts)
{
    std::map<size_t, int> node_rank_map; 
    std::map<int, std::set<size_t>> rank_nodes_map;
    mesh.populate_node_rank_maps(node_rank_map, rank_nodes_map, mesh_maps.first, 2);
    std::map<size_t, std::set<size_t>> node_element_map;
    mesh.populate_node_element_map(node_element_map, mesh_maps.second);

    std::set<size_t> elem_id_set_on_rank_0;
    mesh.find_rank_elements(elem_id_set_on_rank_0, rank_nodes_map, node_element_map, 0);
    int rank0_nelems = elem_id_set_on_rank_0.size();

    std::set<size_t> elem_id_set_on_rank_1;
    mesh.find_rank_elements(elem_id_set_on_rank_1, rank_nodes_map, node_element_map, 1);
    int rank1_nelems = elem_id_set_on_rank_1.size();

    ElemIdNodeIdPairVector elem_nodes_vector_on_rank_0;
    elem_nodes_vector_on_rank_0.reserve(rank0_nelems);
    mesh.filter_element_vector(elem_nodes_vector_on_rank_0, elem_id_set_on_rank_0, mesh_maps.second, 0);

    ElemIdNodeIdPairVector elem_nodes_vector_on_rank_1;
    elem_nodes_vector_on_rank_1.reserve(rank1_nelems);
    mesh.filter_element_vector(elem_nodes_vector_on_rank_1, elem_id_set_on_rank_1, mesh_maps.second, 1);

    std::set<size_t> node_id_set_on_rank_0;
    mesh.find_rank_nodes(node_id_set_on_rank_0, elem_nodes_vector_on_rank_0);
    int rank0_nnodes = node_id_set_on_rank_0.size();

    std::set<size_t> node_id_set_on_rank_1;
    mesh.find_rank_nodes(node_id_set_on_rank_1, elem_nodes_vector_on_rank_1);
    int rank1_nnodes = node_id_set_on_rank_1.size();

    EXPECT_EQ(rank0_nnodes, 6);
    EXPECT_EQ(rank1_nnodes, 6);
}

TEST_F(DistributedLineMeshTests, find_rank_nodes_10_on_2_values)
{
    std::map<size_t, int> node_rank_map; 
    std::map<int, std::set<size_t>> rank_nodes_map;
    mesh.populate_node_rank_maps(node_rank_map, rank_nodes_map, mesh_maps.first, 2);
    std::map<size_t, std::set<size_t>> node_element_map;
    mesh.populate_node_element_map(node_element_map, mesh_maps.second);

    std::set<size_t> elem_id_set_on_rank_0;
    mesh.find_rank_elements(elem_id_set_on_rank_0, rank_nodes_map, node_element_map, 0);
    int rank0_nelems = elem_id_set_on_rank_0.size();

    std::set<size_t> elem_id_set_on_rank_1;
    mesh.find_rank_elements(elem_id_set_on_rank_1, rank_nodes_map, node_element_map, 1);
    int rank1_nelems = elem_id_set_on_rank_1.size();

    ElemIdNodeIdPairVector elem_nodes_vector_on_rank_0;
    elem_nodes_vector_on_rank_0.reserve(rank0_nelems);
    mesh.filter_element_vector(elem_nodes_vector_on_rank_0, elem_id_set_on_rank_0, mesh_maps.second, 0);

    ElemIdNodeIdPairVector elem_nodes_vector_on_rank_1;
    elem_nodes_vector_on_rank_1.reserve(rank1_nelems);
    mesh.filter_element_vector(elem_nodes_vector_on_rank_1, elem_id_set_on_rank_1, mesh_maps.second, 1);

    std::set<size_t> node_id_set_on_rank_0;
    mesh.find_rank_nodes(node_id_set_on_rank_0, elem_nodes_vector_on_rank_0);
    int rank0_nnodes = node_id_set_on_rank_0.size();

    std::set<size_t> node_id_set_on_rank_1;
    mesh.find_rank_nodes(node_id_set_on_rank_1, elem_nodes_vector_on_rank_1);
    int rank1_nnodes = node_id_set_on_rank_1.size();

    check_set_ids(node_id_set_on_rank_0, 0, 6, 1);
    check_set_ids(node_id_set_on_rank_1, 0, 6, 5);
}

TEST_F(DistributedLineMeshTests, find_rank_nodes_10_on_3_counts)
{
    std::map<size_t, int> node_rank_map; 
    std::map<int, std::set<size_t>> rank_nodes_map;
    mesh.populate_node_rank_maps(node_rank_map, rank_nodes_map, mesh_maps.first, 3);
    std::map<size_t, std::set<size_t>> node_element_map;
    mesh.populate_node_element_map(node_element_map, mesh_maps.second);

    std::set<size_t> elem_id_set_on_rank_0;
    mesh.find_rank_elements(elem_id_set_on_rank_0, rank_nodes_map, node_element_map, 0);
    int rank0_nelems = elem_id_set_on_rank_0.size();

    std::set<size_t> elem_id_set_on_rank_1;
    mesh.find_rank_elements(elem_id_set_on_rank_1, rank_nodes_map, node_element_map, 1);
    int rank1_nelems = elem_id_set_on_rank_1.size();

    std::set<size_t> elem_id_set_on_rank_2;
    mesh.find_rank_elements(elem_id_set_on_rank_2, rank_nodes_map, node_element_map, 2);
    int rank2_nelems = elem_id_set_on_rank_2.size();

    ElemIdNodeIdPairVector elem_nodes_vector_on_rank_0;
    elem_nodes_vector_on_rank_0.reserve(rank0_nelems);
    mesh.filter_element_vector(elem_nodes_vector_on_rank_0, elem_id_set_on_rank_0, mesh_maps.second, 0);

    ElemIdNodeIdPairVector elem_nodes_vector_on_rank_1;
    elem_nodes_vector_on_rank_1.reserve(rank1_nelems);
    mesh.filter_element_vector(elem_nodes_vector_on_rank_1, elem_id_set_on_rank_1, mesh_maps.second, 1);

    ElemIdNodeIdPairVector elem_nodes_vector_on_rank_2;
    elem_nodes_vector_on_rank_2.reserve(rank2_nelems);
    mesh.filter_element_vector(elem_nodes_vector_on_rank_2, elem_id_set_on_rank_2, mesh_maps.second, 2);

    std::set<size_t> node_id_set_on_rank_0;
    mesh.find_rank_nodes(node_id_set_on_rank_0, elem_nodes_vector_on_rank_0);
    int rank0_nnodes = node_id_set_on_rank_0.size();

    std::set<size_t> node_id_set_on_rank_1;
    mesh.find_rank_nodes(node_id_set_on_rank_1, elem_nodes_vector_on_rank_1);
    int rank1_nnodes = node_id_set_on_rank_1.size();

    std::set<size_t> node_id_set_on_rank_2;
    mesh.find_rank_nodes(node_id_set_on_rank_2, elem_nodes_vector_on_rank_2);
    int rank2_nnodes = node_id_set_on_rank_2.size();

    EXPECT_EQ(rank0_nnodes, 4);
    EXPECT_EQ(rank1_nnodes, 5);
    EXPECT_EQ(rank2_nnodes, 5);
}

TEST_F(DistributedLineMeshTests, find_rank_nodes_10_on_3_values)
{
    std::map<size_t, int> node_rank_map; 
    std::map<int, std::set<size_t>> rank_nodes_map;
    mesh.populate_node_rank_maps(node_rank_map, rank_nodes_map, mesh_maps.first, 3);
    std::map<size_t, std::set<size_t>> node_element_map;
    mesh.populate_node_element_map(node_element_map, mesh_maps.second);

    std::set<size_t> elem_id_set_on_rank_0;
    mesh.find_rank_elements(elem_id_set_on_rank_0, rank_nodes_map, node_element_map, 0);
    int rank0_nelems = elem_id_set_on_rank_0.size();

    std::set<size_t> elem_id_set_on_rank_1;
    mesh.find_rank_elements(elem_id_set_on_rank_1, rank_nodes_map, node_element_map, 1);
    int rank1_nelems = elem_id_set_on_rank_1.size();

    std::set<size_t> elem_id_set_on_rank_2;
    mesh.find_rank_elements(elem_id_set_on_rank_2, rank_nodes_map, node_element_map, 2);
    int rank2_nelems = elem_id_set_on_rank_2.size();

    ElemIdNodeIdPairVector elem_nodes_vector_on_rank_0;
    elem_nodes_vector_on_rank_0.reserve(rank0_nelems);
    mesh.filter_element_vector(elem_nodes_vector_on_rank_0, elem_id_set_on_rank_0, mesh_maps.second, 0);

    ElemIdNodeIdPairVector elem_nodes_vector_on_rank_1;
    elem_nodes_vector_on_rank_1.reserve(rank1_nelems);
    mesh.filter_element_vector(elem_nodes_vector_on_rank_1, elem_id_set_on_rank_1, mesh_maps.second, 1);

    ElemIdNodeIdPairVector elem_nodes_vector_on_rank_2;
    elem_nodes_vector_on_rank_2.reserve(rank2_nelems);
    mesh.filter_element_vector(elem_nodes_vector_on_rank_2, elem_id_set_on_rank_2, mesh_maps.second, 2);
    std::set<size_t> node_id_set_on_rank_0;
    mesh.find_rank_nodes(node_id_set_on_rank_0, elem_nodes_vector_on_rank_0);
    int rank0_nnodes = node_id_set_on_rank_0.size();

    std::set<size_t> node_id_set_on_rank_1;
    mesh.find_rank_nodes(node_id_set_on_rank_1, elem_nodes_vector_on_rank_1);
    int rank1_nnodes = node_id_set_on_rank_1.size();

    std::set<size_t> node_id_set_on_rank_2;
    mesh.find_rank_nodes(node_id_set_on_rank_2, elem_nodes_vector_on_rank_2);
    int rank2_nnodes = node_id_set_on_rank_2.size();

    check_set_ids(node_id_set_on_rank_0, 0, 4, 1);
    check_set_ids(node_id_set_on_rank_1, 0, 5, 3);
    check_set_ids(node_id_set_on_rank_2, 0, 5, 6);
}
//@}
/**
 * @name filter_node_vector
 * @brief Tests the filter_node_vector function with different num_ranks and for each rank. The IDs should correctly be those in the node ID set provided.
 */
//@{
TEST_F(DistributedLineMeshTests, filter_node_vector_10_on_1_counts)
{
    std::map<size_t, int> node_rank_map; 
    std::map<int, std::set<size_t>> rank_nodes_map;
    mesh.populate_node_rank_maps(node_rank_map, rank_nodes_map, mesh_maps.first, 1);
    std::map<size_t, std::set<size_t>> node_element_map;
    mesh.populate_node_element_map(node_element_map, mesh_maps.second);

    std::set<size_t> elem_id_set_on_rank_0;
    mesh.find_rank_elements(elem_id_set_on_rank_0, rank_nodes_map, node_element_map, 0);
    int rank_nelems = elem_id_set_on_rank_0.size();

    ElemIdNodeIdPairVector elem_nodes_vector_on_rank_0;
    elem_nodes_vector_on_rank_0.reserve(rank_nelems);
    mesh.filter_element_vector(elem_nodes_vector_on_rank_0, elem_id_set_on_rank_0, mesh_maps.second, 0);

    std::set<size_t> node_id_set_on_rank_0;
    mesh.find_rank_nodes(node_id_set_on_rank_0, elem_nodes_vector_on_rank_0);
    int rank0_nnodes = node_id_set_on_rank_0.size();

    NodeIdCoordsPairsVector nodes_coords_vector_on_rank_0;
    nodes_coords_vector_on_rank_0.reserve(rank0_nnodes);
    mesh.filter_node_vector(nodes_coords_vector_on_rank_0, mesh_maps.first, node_id_set_on_rank_0, 0);

    EXPECT_EQ(nodes_coords_vector_on_rank_0.size(), 10);
}

TEST_F(DistributedLineMeshTests, filter_node_vector_10_on_1_values)
{
    std::map<size_t, int> node_rank_map; 
    std::map<int, std::set<size_t>> rank_nodes_map;
    mesh.populate_node_rank_maps(node_rank_map, rank_nodes_map, mesh_maps.first, 1);
    std::map<size_t, std::set<size_t>> node_element_map;
    mesh.populate_node_element_map(node_element_map, mesh_maps.second);

    std::set<size_t> elem_id_set_on_rank_0;
    mesh.find_rank_elements(elem_id_set_on_rank_0, rank_nodes_map, node_element_map, 0);
    int rank0_nelems = elem_id_set_on_rank_0.size();

    ElemIdNodeIdPairVector elem_nodes_vector_on_rank_0;
    elem_nodes_vector_on_rank_0.reserve(rank0_nelems);
    mesh.filter_element_vector(elem_nodes_vector_on_rank_0, elem_id_set_on_rank_0, mesh_maps.second, 0);

    std::set<size_t> node_id_set_on_rank_0;
    mesh.find_rank_nodes(node_id_set_on_rank_0, elem_nodes_vector_on_rank_0);
    int rank0_nnodes = node_id_set_on_rank_0.size();

    NodeIdCoordsPairsVector nodes_coords_vector_on_rank_0;
    nodes_coords_vector_on_rank_0.reserve(rank0_nnodes);
    mesh.filter_node_vector(nodes_coords_vector_on_rank_0, mesh_maps.first, node_id_set_on_rank_0, 0);

    check_vector_ids(nodes_coords_vector_on_rank_0, 0, 10, 1);
}

TEST_F(DistributedLineMeshTests, filter_node_vector_10_on_2_counts)
{
    std::map<size_t, int> node_rank_map; 
    std::map<int, std::set<size_t>> rank_nodes_map;
    mesh.populate_node_rank_maps(node_rank_map, rank_nodes_map, mesh_maps.first, 2);
    std::map<size_t, std::set<size_t>> node_element_map;
    mesh.populate_node_element_map(node_element_map, mesh_maps.second);

    std::set<size_t> elem_id_set_on_rank_0;
    mesh.find_rank_elements(elem_id_set_on_rank_0, rank_nodes_map, node_element_map, 0);
    int rank0_nelems = elem_id_set_on_rank_0.size();

    std::set<size_t> elem_id_set_on_rank_1;
    mesh.find_rank_elements(elem_id_set_on_rank_1, rank_nodes_map, node_element_map, 1);
    int rank1_nelems = elem_id_set_on_rank_1.size();

    ElemIdNodeIdPairVector elem_nodes_vector_on_rank_0;
    elem_nodes_vector_on_rank_0.reserve(rank0_nelems);
    mesh.filter_element_vector(elem_nodes_vector_on_rank_0, elem_id_set_on_rank_0, mesh_maps.second, 0);

    ElemIdNodeIdPairVector elem_nodes_vector_on_rank_1;
    elem_nodes_vector_on_rank_1.reserve(rank1_nelems);
    mesh.filter_element_vector(elem_nodes_vector_on_rank_1, elem_id_set_on_rank_1, mesh_maps.second, 1);

    std::set<size_t> node_id_set_on_rank_0;
    mesh.find_rank_nodes(node_id_set_on_rank_0, elem_nodes_vector_on_rank_0);
    int rank0_nnodes = node_id_set_on_rank_0.size();

    std::set<size_t> node_id_set_on_rank_1;
    mesh.find_rank_nodes(node_id_set_on_rank_1, elem_nodes_vector_on_rank_1);
    int rank1_nnodes = node_id_set_on_rank_1.size();

    NodeIdCoordsPairsVector nodes_coords_vector_on_rank_0;
    nodes_coords_vector_on_rank_0.reserve(rank0_nnodes);
    mesh.filter_node_vector(nodes_coords_vector_on_rank_0, mesh_maps.first, node_id_set_on_rank_0, 0);

    NodeIdCoordsPairsVector nodes_coords_vector_on_rank_1;
    nodes_coords_vector_on_rank_1.reserve(rank1_nnodes);
    mesh.filter_node_vector(nodes_coords_vector_on_rank_1, mesh_maps.first, node_id_set_on_rank_1, 1);

    EXPECT_EQ(nodes_coords_vector_on_rank_0.size(), 6);
    EXPECT_EQ(nodes_coords_vector_on_rank_1.size(), 6);
}

TEST_F(DistributedLineMeshTests, filter_node_vector_10_on_2_values)
{
    std::map<size_t, int> node_rank_map; 
    std::map<int, std::set<size_t>> rank_nodes_map;
    mesh.populate_node_rank_maps(node_rank_map, rank_nodes_map, mesh_maps.first, 2);
    std::map<size_t, std::set<size_t>> node_element_map;
    mesh.populate_node_element_map(node_element_map, mesh_maps.second);

    std::set<size_t> elem_id_set_on_rank_0;
    mesh.find_rank_elements(elem_id_set_on_rank_0, rank_nodes_map, node_element_map, 0);
    int rank0_nelems = elem_id_set_on_rank_0.size();

    std::set<size_t> elem_id_set_on_rank_1;
    mesh.find_rank_elements(elem_id_set_on_rank_1, rank_nodes_map, node_element_map, 1);
    int rank1_nelems = elem_id_set_on_rank_1.size();

    ElemIdNodeIdPairVector elem_nodes_vector_on_rank_0;
    elem_nodes_vector_on_rank_0.reserve(rank0_nelems);
    mesh.filter_element_vector(elem_nodes_vector_on_rank_0, elem_id_set_on_rank_0, mesh_maps.second, 0);

    ElemIdNodeIdPairVector elem_nodes_vector_on_rank_1;
    elem_nodes_vector_on_rank_1.reserve(rank1_nelems);
    mesh.filter_element_vector(elem_nodes_vector_on_rank_1, elem_id_set_on_rank_1, mesh_maps.second, 1);

    std::set<size_t> node_id_set_on_rank_0;
    mesh.find_rank_nodes(node_id_set_on_rank_0, elem_nodes_vector_on_rank_0);
    int rank0_nnodes = node_id_set_on_rank_0.size();

    std::set<size_t> node_id_set_on_rank_1;
    mesh.find_rank_nodes(node_id_set_on_rank_1, elem_nodes_vector_on_rank_1);
    int rank1_nnodes = node_id_set_on_rank_1.size();

    NodeIdCoordsPairsVector nodes_coords_vector_on_rank_0;
    nodes_coords_vector_on_rank_0.reserve(rank0_nnodes);
    mesh.filter_node_vector(nodes_coords_vector_on_rank_0, mesh_maps.first, node_id_set_on_rank_0, 0);

    NodeIdCoordsPairsVector nodes_coords_vector_on_rank_1;
    nodes_coords_vector_on_rank_1.reserve(rank1_nnodes);
    mesh.filter_node_vector(nodes_coords_vector_on_rank_1, mesh_maps.first, node_id_set_on_rank_1, 1);

    check_vector_ids(nodes_coords_vector_on_rank_0, 0, 6, 1);
    check_vector_ids(nodes_coords_vector_on_rank_1, 0, 6, 5);
}

TEST_F(DistributedLineMeshTests, filter_node_vector_10_on_3_counts)
{
    std::map<size_t, int> node_rank_map; 
    std::map<int, std::set<size_t>> rank_nodes_map;
    mesh.populate_node_rank_maps(node_rank_map, rank_nodes_map, mesh_maps.first, 3);
    std::map<size_t, std::set<size_t>> node_element_map;
    mesh.populate_node_element_map(node_element_map, mesh_maps.second);

    std::set<size_t> elem_id_set_on_rank_0;
    mesh.find_rank_elements(elem_id_set_on_rank_0, rank_nodes_map, node_element_map, 0);
    int rank0_nelems = elem_id_set_on_rank_0.size();

    std::set<size_t> elem_id_set_on_rank_1;
    mesh.find_rank_elements(elem_id_set_on_rank_1, rank_nodes_map, node_element_map, 1);
    int rank1_nelems = elem_id_set_on_rank_1.size();

    std::set<size_t> elem_id_set_on_rank_2;
    mesh.find_rank_elements(elem_id_set_on_rank_2, rank_nodes_map, node_element_map, 2);
    int rank2_nelems = elem_id_set_on_rank_2.size();

    ElemIdNodeIdPairVector elem_nodes_vector_on_rank_0;
    elem_nodes_vector_on_rank_0.reserve(rank0_nelems);
    mesh.filter_element_vector(elem_nodes_vector_on_rank_0, elem_id_set_on_rank_0, mesh_maps.second, 0);

    ElemIdNodeIdPairVector elem_nodes_vector_on_rank_1;
    elem_nodes_vector_on_rank_1.reserve(rank1_nelems);
    mesh.filter_element_vector(elem_nodes_vector_on_rank_1, elem_id_set_on_rank_1, mesh_maps.second, 1);

    ElemIdNodeIdPairVector elem_nodes_vector_on_rank_2;
    elem_nodes_vector_on_rank_2.reserve(rank2_nelems);
    mesh.filter_element_vector(elem_nodes_vector_on_rank_2, elem_id_set_on_rank_2, mesh_maps.second, 2);

    std::set<size_t> node_id_set_on_rank_0;
    mesh.find_rank_nodes(node_id_set_on_rank_0, elem_nodes_vector_on_rank_0);
    int rank0_nnodes = node_id_set_on_rank_0.size();

    std::set<size_t> node_id_set_on_rank_1;
    mesh.find_rank_nodes(node_id_set_on_rank_1, elem_nodes_vector_on_rank_1);
    int rank1_nnodes = node_id_set_on_rank_1.size();

    std::set<size_t> node_id_set_on_rank_2;
    mesh.find_rank_nodes(node_id_set_on_rank_2, elem_nodes_vector_on_rank_2);
    int rank2_nnodes = node_id_set_on_rank_2.size();

    NodeIdCoordsPairsVector nodes_coords_vector_on_rank_0;
    nodes_coords_vector_on_rank_0.reserve(rank0_nnodes);
    mesh.filter_node_vector(nodes_coords_vector_on_rank_0, mesh_maps.first, node_id_set_on_rank_0, 0);

    NodeIdCoordsPairsVector nodes_coords_vector_on_rank_1;
    nodes_coords_vector_on_rank_1.reserve(rank1_nnodes);
    mesh.filter_node_vector(nodes_coords_vector_on_rank_1, mesh_maps.first, node_id_set_on_rank_1, 1);

    NodeIdCoordsPairsVector nodes_coords_vector_on_rank_2;
    nodes_coords_vector_on_rank_2.reserve(rank2_nnodes);
    mesh.filter_node_vector(nodes_coords_vector_on_rank_2, mesh_maps.first, node_id_set_on_rank_2, 2);  

    EXPECT_EQ(nodes_coords_vector_on_rank_0.size(), 4);
    EXPECT_EQ(nodes_coords_vector_on_rank_1.size(), 5);
    EXPECT_EQ(nodes_coords_vector_on_rank_2.size(), 5);
}

TEST_F(DistributedLineMeshTests, filter_node_vector_10_on_3_values)
{
    std::map<size_t, int> node_rank_map; 
    std::map<int, std::set<size_t>> rank_nodes_map;
    mesh.populate_node_rank_maps(node_rank_map, rank_nodes_map, mesh_maps.first, 3);
    std::map<size_t, std::set<size_t>> node_element_map;
    mesh.populate_node_element_map(node_element_map, mesh_maps.second);

    std::set<size_t> elem_id_set_on_rank_0;
    mesh.find_rank_elements(elem_id_set_on_rank_0, rank_nodes_map, node_element_map, 0);
    int rank0_nelems = elem_id_set_on_rank_0.size();

    std::set<size_t> elem_id_set_on_rank_1;
    mesh.find_rank_elements(elem_id_set_on_rank_1, rank_nodes_map, node_element_map, 1);
    int rank1_nelems = elem_id_set_on_rank_1.size();

    std::set<size_t> elem_id_set_on_rank_2;
    mesh.find_rank_elements(elem_id_set_on_rank_2, rank_nodes_map, node_element_map, 2);
    int rank2_nelems = elem_id_set_on_rank_2.size();

    ElemIdNodeIdPairVector elem_nodes_vector_on_rank_0;
    elem_nodes_vector_on_rank_0.reserve(rank0_nelems);
    mesh.filter_element_vector(elem_nodes_vector_on_rank_0, elem_id_set_on_rank_0, mesh_maps.second, 0);

    ElemIdNodeIdPairVector elem_nodes_vector_on_rank_1;
    elem_nodes_vector_on_rank_1.reserve(rank1_nelems);
    mesh.filter_element_vector(elem_nodes_vector_on_rank_1, elem_id_set_on_rank_1, mesh_maps.second, 1);

    ElemIdNodeIdPairVector elem_nodes_vector_on_rank_2;
    elem_nodes_vector_on_rank_2.reserve(rank2_nelems);
    mesh.filter_element_vector(elem_nodes_vector_on_rank_2, elem_id_set_on_rank_2, mesh_maps.second, 2);
    std::set<size_t> node_id_set_on_rank_0;
    mesh.find_rank_nodes(node_id_set_on_rank_0, elem_nodes_vector_on_rank_0);
    int rank0_nnodes = node_id_set_on_rank_0.size();

    std::set<size_t> node_id_set_on_rank_1;
    mesh.find_rank_nodes(node_id_set_on_rank_1, elem_nodes_vector_on_rank_1);
    int rank1_nnodes = node_id_set_on_rank_1.size();

    std::set<size_t> node_id_set_on_rank_2;
    mesh.find_rank_nodes(node_id_set_on_rank_2, elem_nodes_vector_on_rank_2);
    int rank2_nnodes = node_id_set_on_rank_2.size();

    NodeIdCoordsPairsVector nodes_coords_vector_on_rank_0;
    nodes_coords_vector_on_rank_0.reserve(rank0_nnodes);
    mesh.filter_node_vector(nodes_coords_vector_on_rank_0, mesh_maps.first, node_id_set_on_rank_0, 0);

    NodeIdCoordsPairsVector nodes_coords_vector_on_rank_1;
    nodes_coords_vector_on_rank_1.reserve(rank1_nnodes);
    mesh.filter_node_vector(nodes_coords_vector_on_rank_1, mesh_maps.first, node_id_set_on_rank_1, 1);

    NodeIdCoordsPairsVector nodes_coords_vector_on_rank_2;
    nodes_coords_vector_on_rank_2.reserve(rank2_nnodes);
    mesh.filter_node_vector(nodes_coords_vector_on_rank_2, mesh_maps.first, node_id_set_on_rank_2, 2);  

    check_vector_ids(nodes_coords_vector_on_rank_0, 0, 4, 1);
    check_vector_ids(nodes_coords_vector_on_rank_1, 0, 5, 3);
    check_vector_ids(nodes_coords_vector_on_rank_2, 0, 5, 6);
}
//@}
#endif 