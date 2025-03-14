#ifndef PLASTIC_MODEL_TESTS_HPP
#define PLASTIC_MODEL_TESTS_HPP

#include "TestHelpers.hpp"

#define PLASTIC_ELEMENT_TYPE NonlinearPlastic


class MeshTestsPlastic : public ::testing::Test {
  public:
    Model model;
    int divisions = 10;
    CommonSectionDefinitions common;
    void SetUp() override {
        common.initialise_section();
        model.create_line_mesh(divisions, {{0.0, 0.0, 0.0}, {10.0, 0.0, 0.0}}, PLASTIC_ELEMENT_TYPE, common.I_section);
    }
    void TearDown() override {
}
};

TEST_F(MeshTestsPlastic, NumOfDivisions)
{
    int num_elems = model.glob_mesh.get_num_elems();
    EXPECT_EQ(num_elems, divisions);
}
TEST_F(MeshTestsPlastic, NumOfNodes)
{
    int num_nodes = model.glob_mesh.get_num_nodes();
    EXPECT_EQ(num_nodes, divisions + 1);
}

#endif 