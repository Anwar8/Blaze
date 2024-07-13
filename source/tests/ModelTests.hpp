#ifndef MODEL_TESTS_HPP
#define MODEL_TESTS_HPP
#include "gtest/gtest.h"
#include "../Model.hpp"
#include <numeric>
#define LOAD_TOLERANCE 1e-6
#define DISP_TOLERANCE 1e-6
#define SOLUTION_TOLERANCE 1e-4

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
    unsigned tracked_node_id = divisions+1;

    void SetUp() override {
        model.create_line_mesh(divisions, {{0.0, 0.0, 0.0}, {10.0, 0.0, 0.0}});
        model.scribe.track_nodes_by_id(std::set<unsigned>{tracked_node_id}, std::set<int>{tracked_dof}, model.glob_mesh);
        
    }
    void TearDown() override {
        
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
    std::shared_ptr<Node> node = model.glob_mesh.get_node_by_id(tracked_node_id);

    node->set_nodal_displacement(tracked_dof, 1.0);
    model.scribe.write_to_records();

    std::vector<Record> record_library = model.scribe.get_record_library();
    Record record = record_library.back();

    std::array<std::vector<real>, 6> recorded_data = record.get_recorded_data();


    EXPECT_NEAR(recorded_data[tracked_dof].back(), 1.0, DISP_TOLERANCE);
}

TEST_F(ScribeTests, CheckTrackedNodeDispTwice)
{
    std::shared_ptr<Node> node = model.glob_mesh.get_node_by_id(tracked_node_id);

    node->set_nodal_displacement(tracked_dof, 1.0);
    model.scribe.write_to_records();
    node->set_nodal_displacement(tracked_dof, 2.0);
    model.scribe.write_to_records();

    std::vector<Record> record_library = model.scribe.get_record_library();
    Record record = record_library.back();

    std::array<std::vector<real>, 6> recorded_data = record.get_recorded_data(); 
    std::vector<real> disp_data = recorded_data[tracked_dof];


    EXPECT_NEAR(disp_data[0], 1.0, DISP_TOLERANCE);
    EXPECT_NEAR(disp_data[1], 2.0, DISP_TOLERANCE);
}



class CantileverBeam : public ::testing::Test {
  public:
    Model model;
    int divisions = 10;
    real y_load = -1e5;
    unsigned tracked_node_id = divisions + 1;
    int tracked_dof = 2;
    real beam_length = 10.0;

    void SetUp() override {
        model.create_line_mesh(divisions, {{0.0, 0.0, 0.0}, {beam_length, 0.0, 0.0}});

        NodalRestraint end_restraint;
        end_restraint.assign_dofs_restraints(std::set<int>{0, 1, 2, 3, 4, 5});
        end_restraint.assign_nodes_by_id(std::set<int>{1}, model.glob_mesh);
        model.restraints.push_back(end_restraint);
        
        NodalRestraint out_of_plane_restraint; 
        out_of_plane_restraint.assign_dofs_restraints(std::set<int>{1, 3, 4});
        out_of_plane_restraint.assign_nodes_by_id(std::set<int>{2, 3, 4, 5, 6, 7, 8, 9, 10, 11}, model.glob_mesh);end_restraint.assign_nodes_by_id(std::set<int>{1}, model.glob_mesh);
        model.restraints.push_back(out_of_plane_restraint);

        model.load_manager.create_a_nodal_load_by_id({(unsigned)(divisions+1)}, std::set<int>{tracked_dof}, std::vector<real>{y_load}, model.glob_mesh);

        model.scribe.track_nodes_by_id(std::set<unsigned>{tracked_node_id}, std::set<int>{tracked_dof}, model.glob_mesh);

        model.initialise_restraints_n_loads();
        model.glob_mesh.check_nodal_loads();

        model.initialise_solution_parameters(1.0, 100, 1e-3, 30);
        model.solve(1);
    }
    void TearDown() override {
}
};


TEST_F(CantileverBeam, CheckResult)
{
    std::shared_ptr<Node> node = model.glob_mesh.get_node_by_id(tracked_node_id);

    std::vector<Record> record_library = model.scribe.get_record_library();
    Record record = record_library.back();

    std::array<std::vector<real>, 6> recorded_data = record.get_recorded_data(); 
    std::vector<real> disp_data = recorded_data[tracked_dof];
    // $\delta = \frac{PL^3}{3EI} = \frac{1e5 (3)^3}{3(2.06e11)(0.0004570000)} = 0.009560026343183701$
    real correct_disp = y_load*std::powf(beam_length, 3)/(3*(2.06e11)*(0.0004570000));
    
    EXPECT_NEAR(disp_data.back(), correct_disp, SOLUTION_TOLERANCE);
}

class SimplySupported : public ::testing::Test {
  public:
    Model model;
    int divisions = 10;
    real y_load = -1e5;
    unsigned mid_node = (divisions/2) + 1;
    
    int tracked_dof = 2;
    real beam_length = 10.0;

    void SetUp() override {
        model.create_line_mesh(divisions, {{0.0, 0.0, 0.0}, {beam_length, 0.0, 0.0}});

        NodalRestraint end_restraint_1;
        end_restraint_1.assign_dofs_restraints(std::set<int>{0, 1, 2, 3, 4}); // restrain x translation, x rotation, y translation, y rotation, and z translation
        end_restraint_1.assign_nodes_by_id(std::set<int>{1}, model.glob_mesh);
        model.restraints.push_back(end_restraint_1);

        NodalRestraint end_restraint_2;
        end_restraint_2.assign_dofs_restraints(std::set<int>{1, 2, 3, 4}); // restrain x rotation, y translation, y rotation, and z translation
        end_restraint_2.assign_nodes_by_id(std::set<int>{divisions + 1}, model.glob_mesh);
        model.restraints.push_back(end_restraint_2);

        NodalRestraint out_of_plane_restraint; 
        out_of_plane_restraint.assign_dofs_restraints(std::set<int>{1, 3, 4}); // restrain x rotation, y rotation, and z translation
        out_of_plane_restraint.assign_nodes_by_id(std::set<int>{2, 3, 4, 5, 6, 7, 8, 9, 10}, model.glob_mesh);
        model.restraints.push_back(out_of_plane_restraint);

        model.load_manager.create_a_nodal_load_by_id({mid_node}, std::set<int>{tracked_dof}, std::vector<real>{y_load}, model.glob_mesh);

        model.scribe.track_nodes_by_id(std::set<unsigned>{mid_node}, std::set<int>{tracked_dof}, model.glob_mesh);

        model.initialise_restraints_n_loads();
        model.glob_mesh.check_nodal_loads();

        model.initialise_solution_parameters(1.0, 100, 1e-3, 30);
        model.solve(1);
    }
    void TearDown() override {
}
};

TEST_F(SimplySupported, CheckResult)
{
    std::shared_ptr<Node> node = model.glob_mesh.get_node_by_id(mid_node);

    std::vector<Record> record_library = model.scribe.get_record_library();
    Record record = record_library.back();

    std::array<std::vector<real>, 6> recorded_data = record.get_recorded_data(); 
    std::vector<real> disp_data = recorded_data[tracked_dof];
    // $\delta = \frac{PL^3}{48EI} $
    real correct_disp = y_load*std::powf(beam_length, 3)/(48*(2.06e11)*(0.0004570000));
    
    EXPECT_NEAR(disp_data.back(), correct_disp, SOLUTION_TOLERANCE);
}


#endif 