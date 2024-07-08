#ifndef MODEL_TESTS_HPP
#define MODEL_TESTS_HPP
#include "gtest/gtest.h"
#include "../Model.hpp"
#include <numeric>
#define LOAD_TOLERANCE 1e-6
#define DISP_TOLERANCE 1e-6

class MeshTests : public ::testing::Test {
  public:
    Model model;
    int divisions = 10;
  
    void SetUp() override {
        model.create_line_mesh(divisions, {{0.0, 0.0, 0.0}, {10.0, 0.0, 0.0}});
    }
    void TearDown() override {
}
};

TEST_F(MeshTests, NumOfDivisions)
{
    int num_elems = model.glob_mesh.get_num_elems();
    EXPECT_EQ(num_elems, divisions);
}
TEST_F(MeshTests, NumOfNodes)
{
    int num_nodes = model.glob_mesh.get_num_nodes();
    EXPECT_EQ(num_nodes, divisions + 1);
}

class RestraintTests : public ::testing::Test {
  public:
    Model model;
    int divisions = 10;
    
    void SetUp() override {
        NodalRestraint end_restraint;
        NodalRestraint out_of_plane_restraint; 
        model.create_line_mesh(divisions, {{0.0, 0.0, 0.0}, {10.0, 0.0, 0.0}});
        end_restraint.assign_dofs_restraints(std::set<int>{0, 1, 2, 3, 4, 5});
        end_restraint.assign_nodes_by_id(std::set<int>{1}, model.glob_mesh);
        model.restraints.push_back(end_restraint);
        out_of_plane_restraint.assign_dofs_restraints(std::set<int>{1, 3, 4});
        out_of_plane_restraint.assign_nodes_by_id(std::set<int>{2, 3, 4, 5, 6, 7, 8, 9, 10, 11}, model.glob_mesh);end_restraint.assign_nodes_by_id(std::set<int>{1}, model.glob_mesh);
    
        model.restraints.push_back(end_restraint);
        model.restraints.push_back(out_of_plane_restraint);


        model.initialise_restraints_n_loads();
    }
    void TearDown() override {
}
};

TEST_F(RestraintTests, ActiveDoFsRestrainedNode)
{
    std::shared_ptr<Node> node = model.glob_mesh.get_node_by_id(1);
    std::set<int> active_dofs = node->get_active_dofs();
    EXPECT_EQ(active_dofs.size(), 0);
}

TEST_F(RestraintTests, RestrainedDoFsRestrainedNode)
{
    std::shared_ptr<Node> node = model.glob_mesh.get_node_by_id(1);
    std::set<int> inactive_dofs = node->get_inactive_dofs();
    EXPECT_EQ(inactive_dofs.size(), 6);
}

TEST_F(RestraintTests, ActiveDoFsFreeNode)
{
    std::shared_ptr<Node> node = model.glob_mesh.get_node_by_id(2);
    std::set<int> active_dofs = node->get_active_dofs();
    EXPECT_EQ(active_dofs.size(), 3);
}

TEST_F(RestraintTests, RestrainedDoFsFreeNode)
{
    std::shared_ptr<Node> node = model.glob_mesh.get_node_by_id(2);
    std::set<int> inactive_dofs = node->get_inactive_dofs();
    EXPECT_EQ(inactive_dofs.size(), 3);
}

class LoadTests : public ::testing::Test {
  public:
    Model model;
    int divisions = 10;
    real y_load = -1e5;

    void SetUp() override {
        model.create_line_mesh(divisions, {{0.0, 0.0, 0.0}, {10.0, 0.0, 0.0}});
        
        model.load_manager.create_a_nodal_load_by_id({(unsigned)(divisions+1)}, std::set<int>{1}, std::vector<real>{y_load}, model.glob_mesh);
        model.initialise_restraints_n_loads();
        
    }
    void TearDown() override {
        model.load_manager.remove_loads();
}
};

TEST_F(LoadTests, UnloadedNodeLoadedDoFs)
{
    std::shared_ptr<Node> node = model.glob_mesh.get_node_by_id(1);
    std::set<int> loaded_dofs = node->get_loaded_dofs();
    EXPECT_EQ(loaded_dofs.size(), 0);
}

TEST_F(LoadTests, LoadedNodeLoadedDoFs)
{
    std::shared_ptr<Node> node = model.glob_mesh.get_node_by_id(divisions+1);
    std::set<int> inactive_dofs = node->get_loaded_dofs();
    EXPECT_EQ(inactive_dofs.size(), 1);
}

TEST_F(LoadTests, LoadedNodeLoadsByDoF)
{
    std::shared_ptr<Node> node = model.glob_mesh.get_node_by_id(divisions+1);
    
    model.load_manager.increment_loads(1.0);
    
    std::array<real, 6> loads = node->get_loads();
    EXPECT_NEAR(loads[1], y_load, LOAD_TOLERANCE);
}

TEST_F(LoadTests, LoadedNodeLoadsByUnloadedDoF)
{
    std::shared_ptr<Node> node = model.glob_mesh.get_node_by_id(divisions+1);
    
    model.load_manager.increment_loads(1.0);
    
    std::array<real, 6> loads = node->get_loads();
    EXPECT_NEAR(loads[0], 0.0, LOAD_TOLERANCE);
}

TEST_F(LoadTests, LoadedNodeTotalLoads)
{
    std::shared_ptr<Node> node = model.glob_mesh.get_node_by_id(divisions+1);
    
    model.load_manager.increment_loads(1.0);
    
    std::array<real, 6> loads = node->get_loads();
    EXPECT_NEAR(std::accumulate(loads.begin(), loads.end(), 0), y_load, LOAD_TOLERANCE);
}

TEST_F(LoadTests, UnloadedNodeTotalLoads)
{
    std::shared_ptr<Node> node = model.glob_mesh.get_node_by_id(3);
    
    model.load_manager.increment_loads(1.0);
    
    std::array<real, 6> loads = node->get_loads();
    EXPECT_NEAR(std::accumulate(loads.begin(), loads.end(), 0), 0.0, LOAD_TOLERANCE);
}

class ScribeTests : public ::testing::Test {
  public:
    Model model;
    int divisions = 10;
    int tracked_dof = 1;

    void SetUp() override {
        model.create_line_mesh(divisions, {{0.0, 0.0, 0.0}, {10.0, 0.0, 0.0}});
        model.scribe.track_nodes_by_id(std::set<unsigned>{(unsigned)(divisions+1)}, std::set<int>{tracked_dof}, model.glob_mesh);
        
    }
    void TearDown() override {
        model.load_manager.remove_loads();
}
};

TEST_F(ScribeTests, CheckNumRecords)
{
    int num_records = model.scribe.get_record_library().size();
    EXPECT_EQ(num_records, 1);
}

TEST_F(ScribeTests, CheckTrackedDofs)
{
    std::vector<Record> records_library = model.scribe.get_record_library();
    std::set<int> dofs = records_library[0].get_tracked_dofs();
    EXPECT_EQ(dofs.size(), 1);
    EXPECT_TRUE(dofs.contains(tracked_dof));
}

TEST_F(ScribeTests, CheckTrackedNodeId)
{
    std::vector<Record> records_library = model.scribe.get_record_library();
    unsigned node_id = records_library[0].get_tracked_node_id();
    EXPECT_EQ(node_id, divisions+1);
}

TEST_F(ScribeTests, CheckTrackedNodeDisp)
{
    Record record = model.scribe.get_record_library()[0];
    std::shared_ptr<Node> node = record.get_tracked_node();
    node->set_nodal_displacement(tracked_dof, 1.0);
    model.scribe.write_to_records();

    auto last_disp = record.get_recorded_data()(0,0);

    EXPECT_NEAR(last_disp, 1.0, DISP_TOLERANCE);
}

TEST_F(ScribeTests, CheckTrackedNodeDispTwice)
{
    Record record = model.scribe.get_record_library()[0];
    std::shared_ptr<Node> node = record.get_tracked_node();
    
    node->set_nodal_displacement(tracked_dof, 1.0);
    model.scribe.write_to_records();

    node->set_nodal_displacement(tracked_dof, 2.0);
    model.scribe.write_to_records();

    auto last_disp = record.get_recorded_data()(1,0);

    EXPECT_NEAR(last_disp, 2.0, DISP_TOLERANCE);
}


// class SolutionTests : public ::testing::Test {
//   public:
//     Model model;
//     int divisions = 10;
//     real y_load = -1e5;

//     void SetUp() override {
//         NodalRestraint end_restraint;
//         NodalRestraint out_of_plane_restraint; 
//         model.create_line_mesh(divisions, {{0.0, 0.0, 0.0}, {10.0, 0.0, 0.0}});
//         end_restraint.assign_dofs_restraints(std::set<int>{0, 1, 2, 3, 4, 5});
//         end_restraint.assign_nodes_by_id(std::set<int>{1}, model.glob_mesh);
//         model.restraints.push_back(end_restraint);
//         out_of_plane_restraint.assign_dofs_restraints(std::set<int>{1, 3, 4});
//         out_of_plane_restraint.assign_nodes_by_id(std::set<int>{2, 3, 4, 5, 6, 7, 8, 9, 10, 11}, model.glob_mesh);end_restraint.assign_nodes_by_id(std::set<int>{1}, model.glob_mesh);
    
//         model.restraints.push_back(end_restraint);
//         model.restraints.push_back(out_of_plane_restraint);

//         model.load_manager.create_a_nodal_load_by_id({(unsigned)(divisions+1)}, std::set<int>{1}, std::vector<real>{y_load}, model.glob_mesh);

//         model.initialise_restraints_n_loads();
//         model.initialise_solution_parameters(1.0, 1, 1e-5, 10);
//     }
//     void TearDown() override {
// }
// };

#endif 