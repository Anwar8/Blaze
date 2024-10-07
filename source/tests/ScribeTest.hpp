#ifndef SCRIBE_TEST_HPP
#define SCRIBE_TEST_HPP

#include "TestHelpers.hpp"

class ScribeOnlyTests : public ::testing::Test {
  public:
    Scribe scribe;
    coords xyz = {0.0, 0.0, 0.0};
    std::shared_ptr<Node> node = std::make_shared<Node>(1, xyz);
    int tracked_dof = 2;


    void SetUp() override {
        scribe.track_nodes_by_ptr({node}, std::set<int>{tracked_dof});
    }
    void TearDown() override {
}
};

TEST_F(ScribeOnlyTests, CheckLibrarySize)
{
    node->set_nodal_displacement(tracked_dof, 1.0);
    scribe.write_to_records();
    

    std::vector<Record> record_library = scribe.get_record_library();
    
    EXPECT_EQ(record_library.size(), 1);
}

TEST_F(ScribeOnlyTests, CheckTrackedNodeDataSizeIsOne)
{
    // node->set_nodal_displacement(tracked_dof, 1.0);
    scribe.write_to_records();
    

    std::vector<Record> record_library = scribe.get_record_library();
    Record record = record_library.back();
    std::array<std::vector<real>, 6> recorded_data = record.get_recorded_data();
    
    EXPECT_EQ(recorded_data[tracked_dof].size(), 1);
}

TEST_F(ScribeOnlyTests, CheckTrackedNodeDataSizeIsTwo)
{
    // node->set_nodal_displacement(tracked_dof, 1.0);
    scribe.write_to_records();
    scribe.write_to_records();
    

    std::vector<Record> record_library = scribe.get_record_library();
    Record record = record_library.back();
    std::array<std::vector<real>, 6> recorded_data = record.get_recorded_data();
    
    EXPECT_EQ(recorded_data[tracked_dof].size(), 2);
}

TEST_F(ScribeOnlyTests, CheckTrackedNodeDataValue)
{
    node->set_nodal_displacement(tracked_dof, 1.0);
    scribe.write_to_records();
    

    std::vector<Record> record_library = scribe.get_record_library();
    Record record = record_library.back();
    std::array<std::vector<real>, 6> recorded_data = record.get_recorded_data();
    std::vector<real> tracked_dof_vector = recorded_data[tracked_dof];
    real last_disp = tracked_dof_vector.back();
    EXPECT_NEAR(last_disp, 1.0, BASIC_TOLERANCE);
}

TEST_F(ScribeOnlyTests, CheckTrackedNodeDataValueTwice)
{
    node->set_nodal_displacement(tracked_dof, 1.0);
    scribe.write_to_records();
    node->set_nodal_displacement(tracked_dof, 2.0);
    scribe.write_to_records();

    std::vector<Record> record_library = scribe.get_record_library();
    Record record = record_library.back();
    std::array<std::vector<real>, 6> recorded_data = record.get_recorded_data();
    std::vector<real> tracked_dof_vector = recorded_data[tracked_dof];
    
    EXPECT_NEAR(tracked_dof_vector[0], 1.0, BASIC_TOLERANCE);
    EXPECT_NEAR(tracked_dof_vector[1], 2.0, BASIC_TOLERANCE);
}

#endif 