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

class CantileverBeamPlastic : public ::testing::Test {
  public:
    Model model;
    CommonSectionDefinitions common;

    real beam_length = 10.0;
    real max_stress = YIELD_STRENGTH/2;
    real correct_moment = common.moment_of_inertia * max_stress/(common.h/2);
    real y_load = -correct_moment/beam_length;
    
    int divisions = 3;
    
    unsigned tracked_node_id = divisions + 1;
    int tracked_dof = 2;
    


    void SetUp() override {
        common.initialise_section();
        model.create_line_mesh(divisions, {{0.0, 0.0, 0.0}, {beam_length, 0.0, 0.0}}, PLASTIC_ELEMENT_TYPE, common.I_section);

        NodalRestraint end_restraint;
        end_restraint.assign_dofs_restraints(std::set<int>{0, 1, 2, 3, 4, 5});
        end_restraint.assign_nodes_by_id(std::set<int>{1}, model.glob_mesh);
        model.restraints.push_back(end_restraint);
        
        NodalRestraint out_of_plane_restraint; 
        out_of_plane_restraint.assign_dofs_restraints(std::set<int>{1, 3, 4});
        // out_of_plane_restraint.assign_nodes_by_id(std::set<int>{2, 3, 4, 5, 6, 7, 8, 9, 10, 11}, model.glob_mesh);end_restraint.assign_nodes_by_id(std::set<int>{1}, model.glob_mesh);
        out_of_plane_restraint.assign_nodes_by_id(std::set<int>{2, 3, 4}, model.glob_mesh);end_restraint.assign_nodes_by_id(std::set<int>{1}, model.glob_mesh);
        model.restraints.push_back(out_of_plane_restraint);

        model.load_manager.create_a_nodal_load_by_id({(unsigned)(divisions+1)}, std::set<int>{tracked_dof}, std::vector<real>{y_load}, model.glob_mesh);

        model.scribe.track_nodes_by_id(std::set<unsigned>{tracked_node_id}, std::set<int>{tracked_dof}, model.glob_mesh);

        model.initialise_restraints_n_loads();
        model.glob_mesh.check_nodal_loads();

        model.initialise_solution_parameters(1.0, 100, 1e-4, 30);
        model.solve(-1);
    }
    void TearDown() override {
}
};


TEST_F(CantileverBeamPlastic, CheckResult)
{
    std::shared_ptr<Node> node = model.glob_mesh.get_node_by_id(tracked_node_id);

    std::vector<Record> record_library = model.scribe.get_record_library();
    Record record = record_library.back();

    std::array<std::vector<real>, 6> recorded_data = record.get_recorded_data(); 
    std::vector<real> disp_data = recorded_data[tracked_dof];
    // $\delta = \frac{PL^3}{3EI} = \frac{1e5 (3)^3}{3(2.06e11)(0.0004570000)} = 0.009560026343183701$
    real correct_disp = y_load*std::pow(beam_length, 3)/(3*(YOUNGS_MODULUS)*(common.moment_of_inertia));
    real tolerance = std::abs(PERCENT_TOLERANCE*correct_disp);
    EXPECT_NEAR(disp_data.back(), correct_disp, tolerance);
}

class SimplySupportedPlastic : public ::testing::Test {
  public:
    Model model;
    CommonSectionDefinitions common;

    int divisions = 10;
    real y_load = -1e5;
    unsigned mid_node = (divisions/2) + 1;
    
    int tracked_dof = 2;
    real beam_length = 10.0;

    void SetUp() override {
        common.initialise_section();
        model.create_line_mesh(divisions, {{0.0, 0.0, 0.0}, {beam_length, 0.0, 0.0}}, PLASTIC_ELEMENT_TYPE, common.I_section);

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
        model.solve(-1);
    }
    void TearDown() override {
}
};

TEST_F(SimplySupportedPlastic, CheckResult)
{
    std::shared_ptr<Node> node = model.glob_mesh.get_node_by_id(mid_node);

    std::vector<Record> record_library = model.scribe.get_record_library();
    Record record = record_library.back();

    std::array<std::vector<real>, 6> recorded_data = record.get_recorded_data(); 
    std::vector<real> disp_data = recorded_data[tracked_dof];
    // $\delta = \frac{PL^3}{48EI} $
    real correct_disp = y_load*std::pow(beam_length, 3)/(48*(YOUNGS_MODULUS)*(common.moment_of_inertia));
    real tolerance = std::abs(PERCENT_TOLERANCE*correct_disp);
    EXPECT_NEAR(disp_data.back(), correct_disp, tolerance);
}

class SimplySupportedUdlPlastic : public ::testing::Test {
  public:
    Model model;
    CommonSectionDefinitions common;
    int divisions = 100;
    real y_udl = -1e4; // N/m
    unsigned mid_node = (divisions/2) + 1;
    std::vector<unsigned> loaded_nodes = std::vector<unsigned>(divisions - 1);
    int tracked_dof = 2;
    real beam_length = 5.0;

    void SetUp() override {
        common.initialise_section();
        model.create_line_mesh(divisions, {{0.0, 0.0, 0.0}, {beam_length, 0.0, 0.0}}, PLASTIC_ELEMENT_TYPE, common.I_section);

        NodalRestraint end_restraint_1;
        end_restraint_1.assign_dofs_restraints(std::set<int>{0, 1, 2, 3, 4}); // restrain x translation, x rotation, y translation, y rotation, and z translation
        end_restraint_1.assign_nodes_by_id(std::set<int>{1}, model.glob_mesh);
        model.restraints.push_back(end_restraint_1);

        NodalRestraint end_restraint_2;
        end_restraint_2.assign_dofs_restraints(std::set<int>{1, 2, 3, 4}); // restrain x rotation, y translation, y rotation, and z translation
        end_restraint_2.assign_nodes_by_id(std::set<int>{divisions + 1}, model.glob_mesh);
        model.restraints.push_back(end_restraint_2);
        // create the loaded and restrained intermediate nodes
        std::iota(loaded_nodes.begin(), loaded_nodes.end(), 2);

        NodalRestraint out_of_plane_restraint; 
        out_of_plane_restraint.assign_dofs_restraints(std::set<int>{1, 3, 4}); // restrain x rotation, y rotation, and z translation
        out_of_plane_restraint.assign_nodes_by_id(loaded_nodes, model.glob_mesh);
        model.restraints.push_back(out_of_plane_restraint);

        // calculate load
        //
        real y_load = y_udl*beam_length/(divisions - 1);
        model.load_manager.create_a_nodal_load_by_id(loaded_nodes, std::set<int>{tracked_dof}, std::vector<real>{y_load}, model.glob_mesh);

        model.scribe.track_nodes_by_id(std::set<unsigned>{mid_node}, std::set<int>{tracked_dof}, model.glob_mesh);

        model.initialise_restraints_n_loads();
        model.glob_mesh.check_nodal_loads();

        model.initialise_solution_parameters(1.0, 100, 1e-3, 10);
        model.solve(-1);
    }
    void TearDown() override {
}
};

TEST_F(SimplySupportedUdlPlastic, CheckResult)
{
    std::shared_ptr<Node> node = model.glob_mesh.get_node_by_id(mid_node);

    std::vector<Record> record_library = model.scribe.get_record_library();
    Record record = record_library.back();

    std::array<std::vector<real>, 6> recorded_data = record.get_recorded_data(); 
    std::vector<real> disp_data = recorded_data[tracked_dof];
    // $\delta = \frac{5 w L^4}{384 EI} $
    real correct_disp = 5*y_udl*std::pow(beam_length, 4)/(384*(YOUNGS_MODULUS)*(common.moment_of_inertia));
    real tolerance = std::abs(PERCENT_TOLERANCE*correct_disp);
    EXPECT_NEAR(disp_data.back(), correct_disp, tolerance);
}

#endif 