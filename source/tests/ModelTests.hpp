#ifndef MODEL_TESTS_HPP
#define MODEL_TESTS_HPP

#include "TestHelpers.hpp"

#define ELEMENT_TYPE NonlinearElastic

class MeshTests : public ::testing::Test {
  public:
    Model model;
    int divisions = 10;
  
    void SetUp() override {
        BasicSection sect(2.06e11, 0.0125, 0.0004570000);
        model.create_line_mesh(divisions, {{0.0, 0.0, 0.0}, {10.0, 0.0, 0.0}}, ELEMENT_TYPE, sect);
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
        BasicSection sect(2.06e11, 0.0125, 0.0004570000);
        model.create_line_mesh(divisions, {{0.0, 0.0, 0.0}, {10.0, 0.0, 0.0}}, ELEMENT_TYPE, sect);
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
        BasicSection sect(2.06e11, 0.0125, 0.0004570000);
        model.create_line_mesh(divisions, {{0.0, 0.0, 0.0}, {10.0, 0.0, 0.0}}, ELEMENT_TYPE, sect);
        
        model.load_manager.create_a_nodal_load_by_id(std::vector<unsigned>{(unsigned)(divisions+1)}, std::set<int>{1}, std::vector<real>{y_load}, model.glob_mesh);
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
    EXPECT_NEAR(loads[1], y_load, BASIC_TOLERANCE);
}

TEST_F(LoadTests, LoadedNodeLoadsByUnloadedDoF)
{
    std::shared_ptr<Node> node = model.glob_mesh.get_node_by_id(divisions+1);
    
    model.load_manager.increment_loads(1.0);
    
    std::array<real, 6> loads = node->get_loads();
    EXPECT_NEAR(loads[0], 0.0, BASIC_TOLERANCE);
}

TEST_F(LoadTests, LoadedNodeTotalLoads)
{
    std::shared_ptr<Node> node = model.glob_mesh.get_node_by_id(divisions+1);
    
    model.load_manager.increment_loads(1.0);
    
    std::array<real, 6> loads = node->get_loads();
    EXPECT_NEAR(std::accumulate(loads.begin(), loads.end(), 0), y_load, BASIC_TOLERANCE);
}

TEST_F(LoadTests, UnloadedNodeTotalLoads)
{
    std::shared_ptr<Node> node = model.glob_mesh.get_node_by_id(3);
    
    model.load_manager.increment_loads(1.0);
    
    std::array<real, 6> loads = node->get_loads();
    EXPECT_NEAR(std::accumulate(loads.begin(), loads.end(), 0), 0.0, BASIC_TOLERANCE);
}

class ScribeTests : public ::testing::Test {
  public:
    Model model;
    int divisions = 10;
    int tracked_dof = 1;
    unsigned tracked_node_id = divisions+1;

    void SetUp() override {
        BasicSection sect(2.06e11, 0.0125, 0.0004570000);
        model.create_line_mesh(divisions, {{0.0, 0.0, 0.0}, {10.0, 0.0, 0.0}}, ELEMENT_TYPE, sect);
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
    EXPECT_TRUE(dofs.count(tracked_dof));
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


    EXPECT_NEAR(recorded_data[tracked_dof].back(), 1.0, BASIC_TOLERANCE);
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


    EXPECT_NEAR(disp_data[0], 1.0, BASIC_TOLERANCE);
    EXPECT_NEAR(disp_data[1], 2.0, BASIC_TOLERANCE);
}




// TEST_F(MacNealSlenderBeam, CheckResultK75)
// {
//     int steps = 75;
//     load_and_run(model, tracked_node_id, loaded_dof, M_max*0.75, steps);
//     std::shared_ptr<Node> node = model.glob_mesh.get_node_by_id(tracked_node_id);

//     std::vector<Record> record_library = model.scribe.get_record_library();
//     Record record = record_library.back();

//     std::array<std::vector<real>, 6> recorded_data = record.get_recorded_data(); 
//     std::vector<real> disp_data = recorded_data[tracked_dof];
//     real disp = disp_data.back();
//     real min_correct_disp = -2.0; // beam should exceed this displacement
//     real max_correct_disp = -2.3; // beam should not exceed this displacement
//     EXPECT_TRUE((-disp) > (-min_correct_disp));
//     EXPECT_TRUE((-disp) < (-max_correct_disp));
//     model.read_all_records();
// }

// TEST_F(MacNealSlenderBeam, CheckResultK100)
// {
//     int steps = 100;
//     load_and_run(model, tracked_node_id, loaded_dof, M_max*1.00, steps);
//     std::shared_ptr<Node> node = model.glob_mesh.get_node_by_id(tracked_node_id);

//     std::vector<Record> record_library = model.scribe.get_record_library();
//     Record record = record_library.back();

//     std::array<std::vector<real>, 6> recorded_data = record.get_recorded_data(); 
//     std::vector<real> disp_data = recorded_data[tracked_dof];
//     real disp = disp_data.back();
//     real min_correct_disp = -0.05; // beam should exceed this displacement
//     real max_correct_disp = 0.05; // beam should not exceed this displacement
//     EXPECT_TRUE((-disp) > (-min_correct_disp));
//     EXPECT_TRUE((-disp) < (-max_correct_disp));
//     model.read_all_records();
// }

#endif 
