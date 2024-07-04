#include "gtest/gtest.h"
#include "../NodalLoad.hpp"
#define NODAL_LOAD_TOL 1e-6

void setup_nodal_load(NodalLoad& nodal_load, std::vector<std::shared_ptr<Node>> in_nodes, std::vector<int> dofs, std::vector<real> loads)
{
    nodal_load.assign_nodes_by_ptr(in_nodes);
    nodal_load.assign_dofs_loads(dofs, loads);
}
void setup_nodal_load(NodalLoad& nodal_load, std::vector<std::shared_ptr<Node>> in_nodes, std::set<int> dofs, std::vector<real> loads)
{
    nodal_load.assign_nodes_by_ptr(in_nodes);
    nodal_load.assign_dofs_loads(dofs, loads);
}

class NodalLoadTest : public ::testing::Test {
    // Declare variables to be used in the fixture
  public:
    std::vector<std::shared_ptr<Node>> in_nodes;
    NodalLoad nodal_load;    
  
    void SetUp() override {
        // Create a set of nodes to be loaded.
        in_nodes.push_back(std::make_shared<Node>(0.0, 0.0, 0.0));
        in_nodes.push_back(std::make_shared<Node>(1.0, 0.0, 0.0));
        in_nodes.push_back(std::make_shared<Node>(2.0, 0.0, 0.0));
    }
    void TearDown() override {
        nodal_load.unload_loaded_nodes();
        nodal_load.reset();
}
};

TEST_F(NodalLoadTest, DofVectorSize)
{
    std::vector<int> dofs = {0, 1, 2};
    std::vector<real> loads = {1.0, 2.0, 3.0};
    setup_nodal_load(nodal_load, in_nodes, dofs, loads);

    EXPECT_EQ(nodal_load.get_loaded_dofs().size(), 3);
}

TEST_F(NodalLoadTest, DofVectorSizeZero)
{
    std::vector<int> dofs = {};
    std::vector<real> loads = {};
    setup_nodal_load(nodal_load, in_nodes, dofs, loads);

    EXPECT_EQ(nodal_load.get_loaded_dofs().size(), 0);
}

TEST_F(NodalLoadTest, ConstructBySetDofSize)
{
    std::set<int> dofs = {0, 1, 2};
    std::vector<real> loads = {1.0, 2.0, 3.0};
    setup_nodal_load(nodal_load, in_nodes, dofs, loads);

    EXPECT_EQ(nodal_load.get_loaded_dofs().size(), 3);
}

TEST_F(NodalLoadTest, LoadOrderCorrect)
{
    std::vector<int> dofs = {0, 5, 2, 3};
    std::vector<real> loads = {1.0, 23.0, -3.0, 1e4};
    setup_nodal_load(nodal_load, in_nodes, dofs, loads);
    auto load_object_loads = nodal_load.get_nodal_loads();

    int i = 0;
    for (auto dof: dofs)
    {
         EXPECT_NEAR(load_object_loads[dof], loads[i], NODAL_LOAD_TOL);
         ++i;
    }
}

TEST_F(NodalLoadTest, ConstructBySetLoadOrderCorrect)
{
    std::set<int> dofs = {0, 5, 2, 3};
    std::vector<real> loads = {1.0, 23.0, -3.0, 1e4};
    setup_nodal_load(nodal_load, in_nodes, dofs, loads);
    auto load_object_loads = nodal_load.get_nodal_loads();

    int i = 0;
    for (auto dof: dofs)
    {
         EXPECT_NEAR(load_object_loads[dof], loads[i], NODAL_LOAD_TOL);
         ++i;
    }
}