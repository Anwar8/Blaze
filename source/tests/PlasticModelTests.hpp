#ifndef PLASTIC_MODEL_TESTS_HPP
#define PLASTIC_MODEL_TESTS_HPP
#include "gtest/gtest.h"
#include "../Model.hpp"
#include "../main.hpp"
#include "ModelTests.hpp"
#include "PlasticBeamElementTests.hpp"
#include <numeric>
#define LOAD_TOLERANCE 1e-6
#define DISP_TOLERANCE 1e-6
#define SOLUTION_TOLERANCE_PERCENT 0.02
#define PLASTIC_ELEMENT_TYPE NonlinearPlastic

struct CommonSectionDefinitions {
    ElasticPlasticMaterial steel = ElasticPlasticMaterial(PLASTIC_YOUNGS_MODULUS_BEAM_TEST, PLASTIC_YIELD_STRENGTH_BEAM_TEST, 0.0);
    BeamColumnFiberSection I_section;
    real tf = 19.6e-3;
    real tw = 11.4e-3;
    real b = 192.8e-3;
    real h = 467.2e-3;
    real d = h - 2*tf;
    int flange_divisions = 10;
    int web_divisions = 40;
    real correct_area = tf*b*2 + (h - 2*tf)*tw; // m2
    real moment_of_inertia = tw*pow(h - 2*tf, 3)/12 + 2*b*pow(tf,3)/12 + 2*(tf*b)*pow(0.5*h - 0.5*tf, 2);
    real section_modulus = moment_of_inertia/(h/2);
    real correct_elastic_moment = section_modulus * PLASTIC_YIELD_STRENGTH_BEAM_TEST;
    real correct_plastic_moment = PLASTIC_YIELD_STRENGTH_BEAM_TEST*(tf*b)*(h - tf) + PLASTIC_YIELD_STRENGTH_BEAM_TEST*((0.5*h - tf)*tw)*(0.5*d);
    real kappa_elastic = correct_elastic_moment/(PLASTIC_YOUNGS_MODULUS_BEAM_TEST * moment_of_inertia);
    real distance_to_first_fibre = (d/40)*0.5;
    real kappa_plastic = PLASTIC_YIELD_STRENGTH_BEAM_TEST/(PLASTIC_YOUNGS_MODULUS_BEAM_TEST*distance_to_first_fibre);
};

class MeshTestsPlastic : public ::testing::Test {
  public:
    Model model;
    int divisions = 10;
    CommonSectionDefinitions common;
    void SetUp() override {
        build_plastic_I_section(common.I_section, common.steel, 0.0, common.tf, common.b, common.tw, common.h, common.flange_divisions, common.web_divisions);
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
    int divisions = 10;
    real y_load = -1e5;
    unsigned tracked_node_id = divisions + 1;
    int tracked_dof = 2;
    real beam_length = 10.0;

    CommonSectionDefinitions common;

    void SetUp() override {
        build_plastic_I_section(common.I_section, common.steel, 0.0, common.tf, common.b, common.tw, common.h, common.flange_divisions, common.web_divisions);
        model.create_line_mesh(divisions, {{0.0, 0.0, 0.0}, {beam_length, 0.0, 0.0}}, LinearElastic, common.I_section);

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
    real correct_disp = y_load*std::powf(beam_length, 3)/(3*(PLASTIC_YOUNGS_MODULUS_BEAM_TEST)*(common.moment_of_inertia));
    real tolerance = std::abs(SOLUTION_TOLERANCE_PERCENT*correct_disp);
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
        build_plastic_I_section(common.I_section, common.steel, 0.0, common.tf, common.b, common.tw, common.h, common.flange_divisions, common.web_divisions);
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
    real correct_disp = y_load*std::powf(beam_length, 3)/(48*(PLASTIC_YOUNGS_MODULUS_BEAM_TEST)*(common.moment_of_inertia));
    real tolerance = std::abs(SOLUTION_TOLERANCE_PERCENT*correct_disp);
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
        build_plastic_I_section(common.I_section, common.steel, 0.0, common.tf, common.b, common.tw, common.h, common.flange_divisions, common.web_divisions);
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

        model.initialise_solution_parameters(1.0, 100, 1e-3, 30);
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
    real correct_disp = 5*y_udl*std::powf(beam_length, 4)/(384*(PLASTIC_YOUNGS_MODULUS_BEAM_TEST)*(common.moment_of_inertia));
    real tolerance = std::abs(SOLUTION_TOLERANCE_PERCENT*correct_disp);
    EXPECT_NEAR(disp_data.back(), correct_disp, tolerance);
}

// class MacNealSlenderBeam : public ::testing::Test {
//   public:
//     Model model;
//     int divisions = 200;
    

//     unsigned tracked_node_id = divisions + 1;
//     int loaded_dof = 5;
//     int tracked_dof = 2;
//     real beam_length = 10.0;
//     real youngs_modulus = 1200000.0;
//     real b = 1.0;
//     real h = 0.1;
//     real moment_of_inertia = (b*h*h*h)/12.0;
//     real area = b*h;

//     std::vector<unsigned> restrained_nodes = std::vector<unsigned>(divisions);

//     real M_max = -20*PI;

//     void SetUp() override {
//         BasicSection sect(youngs_modulus, area, moment_of_inertia);
//         model.create_line_mesh(divisions, {{0.0, 0.0, 0.0}, {beam_length, 0.0, 0.0}}, PLASTIC_ELEMENT_TYPE, sect);
        
//         std::iota(restrained_nodes.begin(), restrained_nodes.end(), 2);

//         NodalRestraint end_restraint;
//         end_restraint.assign_dofs_restraints(std::set<int>{0, 1, 2, 3, 4, 5});
//         end_restraint.assign_nodes_by_id(std::set<int>{1}, model.glob_mesh);
//         model.restraints.push_back(end_restraint);
        
//         NodalRestraint out_of_plane_restraint; 
//         out_of_plane_restraint.assign_dofs_restraints(std::set<int>{1, 3, 4});
//         out_of_plane_restraint.assign_nodes_by_id(restrained_nodes, model.glob_mesh);
//         model.restraints.push_back(out_of_plane_restraint);

        

//         model.scribe.track_nodes_by_id(std::set<unsigned>{tracked_node_id}, std::set<int>{tracked_dof}, model.glob_mesh);
        
//     }
//     void TearDown() override {
// }

// };

// void load_and_run(Model& model, unsigned loaded_node, int loaded_dof, real load, int steps)
// {
//     model.load_manager.create_a_nodal_load_by_id({loaded_node}, std::set<int>{loaded_dof}, std::vector<real>{load}, model.glob_mesh);
//     model.initialise_restraints_n_loads();
//     model.glob_mesh.check_nodal_loads();

//     model.initialise_solution_parameters(1.0, steps, 1e-4, 30);
//     model.solve(100);
// }


// TEST_F(MacNealSlenderBeam, CheckResultK25)
// {
//     int steps = 25;
//     load_and_run(model, tracked_node_id, loaded_dof, M_max*0.25, steps);
//     std::shared_ptr<Node> node = model.glob_mesh.get_node_by_id(tracked_node_id);

//     std::vector<Record> record_library = model.scribe.get_record_library();
//     Record record = record_library.back();

//     std::array<std::vector<real>, 6> recorded_data = record.get_recorded_data(); 
//     std::vector<real> disp_data = recorded_data[tracked_dof];
//     real disp = disp_data.back();
//     real min_correct_disp = -6.0; // beam should exceed this displacement
//     real max_correct_disp = -7.0; // beam should not exceed this displacement
//     EXPECT_TRUE((-disp) > (-min_correct_disp));
//     EXPECT_TRUE((-disp) < (-max_correct_disp));
// }

// TEST_F(MacNealSlenderBeam, CheckResultK50)
// {
//     int steps = 50;
//     load_and_run(model, tracked_node_id, loaded_dof, M_max*0.5, steps);
//     std::shared_ptr<Node> node = model.glob_mesh.get_node_by_id(tracked_node_id);

//     std::vector<Record> record_library = model.scribe.get_record_library();
//     Record record = record_library.back();

//     std::array<std::vector<real>, 6> recorded_data = record.get_recorded_data(); 
//     std::vector<real> disp_data = recorded_data[tracked_dof];
//     real disp = disp_data.back();
//     real min_correct_disp = -6.0; // beam should exceed this displacement
//     real max_correct_disp = -7.0; // beam should not exceed this displacement
//     EXPECT_TRUE((-disp) > (-min_correct_disp));
//     EXPECT_TRUE((-disp) < (-max_correct_disp));
// }

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