/**
 * @file PlasticBeamElementTests.hpp
 * @brief tests the nonlinear plastic beam-column element.
 * 
 */
#ifndef PLASTIC_BEAM_ELEMENT_TESTS_HPP
#define PLASTIC_BEAM_ELEMENT_TESTS_HPP

#define PLASTIC_BEAM_TOLERANCE 1e-6
#include <iostream>
#include <vector>
#include "gtest/gtest.h"
#include "../beam_element.hpp"
#include "../BeamElementBaseClass.hpp"
#include "../BeamElementCommonInterface.hpp"
#include "../Linear2DBeamElement.hpp"
#include "../Nonlinear2DBeamElement.hpp"
#include "../Nonlinear2DPlasticBeamElement.hpp"
#include "../maths_defaults.hpp"
#include "../node.hpp"
#include "BeamElementTests.hpp"

#define PLASTIC_BEAM_LENGTH 3.0
#define PLASTIC_YOUNGS_MODULUS_BEAM_TEST 2.06e11
#define PLASTIC_YIELD_STRENGTH_BEAM_TEST 550.0e6
#define A 0.0125
#define I 0.0004570000


void build_plastic_I_section(BeamColumnFiberSection& section, ElasticPlasticMaterial& steel, real offset, real tf, real b, real tw, real h, int flange_divisions, int web_divisions)
{
    std::vector<real> areas;
    std::vector<real> ys;
    
    real y = offset;
    real area = 0.0;
    // bottom flange
    area = b*tf/(flange_divisions);
    y -= 0.5*(tf/flange_divisions);
    for (int i = 0; i < flange_divisions; ++i)
    {
        y += (tf/flange_divisions);
        ys.push_back(y);
        areas.push_back(area);
    }
    // web
    area = (h - 2*tf)*tw / web_divisions;
    y = offset + tf - 0.5*((h - 2*tf)/web_divisions);
    for (int i = flange_divisions; i < flange_divisions + web_divisions; ++i)
    {
        y += ((h - 2*tf)/web_divisions);
        ys.push_back(y);
        areas.push_back(area);
    }
    // top flange
    area = b*tf/(flange_divisions);
    y = offset + (h - tf)  - 0.5*(tf/flange_divisions);
    for (int i = flange_divisions + web_divisions; i < flange_divisions + web_divisions + flange_divisions; ++i)
    {
        y += (tf/flange_divisions);
        ys.push_back(y);
        areas.push_back(area);
    }

    section.add_fibres(&steel, areas, ys);
}

void common_plastic_beam_setup(std::vector<std::shared_ptr<Node>>& in_nodes, BeamColumnFiberSection& sect, std::shared_ptr<Nonlinear2DPlasticBeamElement>& my_beam, vec& U) 
{
    // Create the nodes
    in_nodes.push_back(std::make_shared<Node>(0.0, 0.0, 0.0));
    in_nodes.push_back(std::make_shared<Node>(PLASTIC_BEAM_LENGTH, 0.0, 0.0));
    my_beam = std::make_shared<Nonlinear2DPlasticBeamElement>(0, in_nodes, sect);
    // Create the d vector
    U = make_xd_vec(12);
}

// void move_right(std::vector<std::shared_ptr<Node>>& in_nodes) 
// {
//     in_nodes[0]->set_nodal_displacement(0, 1.0); // node 1 U1
//     in_nodes[1]->set_nodal_displacement(0, 1.0); // node 2 U1
// }
// void move_up(std::vector<std::shared_ptr<Node>>& in_nodes) 
// {
//     in_nodes[0]->set_nodal_displacement(2, 1.0); // node 1 U2
//     in_nodes[1]->set_nodal_displacement(2, 1.0); // node 2 U2
// }
std::pair<real, real> rotate_ccw(std::vector<std::shared_ptr<Node>>& in_nodes, real theta) 
{
    real delta_x = PLASTIC_BEAM_LENGTH/2.0 - std::cos(theta)*(PLASTIC_BEAM_LENGTH/2.0);
    real delta_y = std::sin(theta)*(PLASTIC_BEAM_LENGTH/2.0);

    in_nodes[0]->set_nodal_displacement(0, delta_x); // node 1 U1
    in_nodes[0]->set_nodal_displacement(2, -delta_y); // node 1 U2
    in_nodes[0]->set_nodal_displacement(5, theta); // node 1 U33
    in_nodes[1]->set_nodal_displacement(0, -delta_x); // node 2 U2
    in_nodes[1]->set_nodal_displacement(2, delta_y); // node 2 U2
    in_nodes[1]->set_nodal_displacement(5, theta); // node 2 U33
    return std::pair<real, real> (delta_x, delta_y);
}

void constant_compression(std::vector<std::shared_ptr<Node>>& in_nodes, real delta) 
{
    in_nodes[0]->set_nodal_displacement(0, delta/2); // node 1 U1
    in_nodes[1]->set_nodal_displacement(0, -delta/2); // node 2 U1
}
void constant_tension(std::vector<std::shared_ptr<Node>>& in_nodes, real delta) 
{
    in_nodes[0]->set_nodal_displacement(0, -delta/2); // node 1 U1
    in_nodes[1]->set_nodal_displacement(0, delta/2); // node 2 U1
}
void constant_positive_bending(std::vector<std::shared_ptr<Node>>& in_nodes, real theta) 
{
    in_nodes[0]->set_nodal_displacement(5, -theta); // node 2 U33
    in_nodes[1]->set_nodal_displacement(5, theta); // node 2 U33
}


class PlasticBeamTests : public ::testing::Test {
  public:
    ElasticPlasticMaterial steel = ElasticPlasticMaterial(PLASTIC_YOUNGS_MODULUS_BEAM_TEST, PLASTIC_YIELD_STRENGTH_BEAM_TEST, 0.0);
    BeamColumnFiberSection I_section;
    std::vector<std::shared_ptr<Node>> in_nodes;
    std::shared_ptr<Nonlinear2DPlasticBeamElement> my_beam;
    vec U;
    
    real tf = 19.6e-3;
    real tw = 11.4e-3;
    real b = 192.8e-3;
    real h = 467.2e-3;
    real d = h - 2*tf;
    real correct_area = tf*b*2 + (h - 2*tf)*tw; // m2

    real moment_of_inertia = tw*pow(h - 2*tf, 3)/12 + 2*b*pow(tf,3)/12 + 2*(tf*b)*pow(0.5*h - 0.5*tf, 2);
    real section_modulus = moment_of_inertia/(h/2);
    real correct_elastic_moment = section_modulus * PLASTIC_YIELD_STRENGTH_BEAM_TEST;
    real correct_plastic_moment = PLASTIC_YIELD_STRENGTH_BEAM_TEST*(tf*b)*(h - tf) + PLASTIC_YIELD_STRENGTH_BEAM_TEST*((0.5*h - tf)*tw)*(0.5*d);
    real kappa_elastic = correct_elastic_moment/(PLASTIC_YOUNGS_MODULUS_BEAM_TEST * moment_of_inertia);

    real distance_to_first_fibre = (d/40)*0.5;
    real kappa_plastic = PLASTIC_YIELD_STRENGTH_BEAM_TEST/(PLASTIC_YOUNGS_MODULUS_BEAM_TEST*distance_to_first_fibre);
    void SetUp() override {
        build_plastic_I_section(I_section, steel, 0.0, tf, b, tw, h, 10, 40);
        common_plastic_beam_setup(in_nodes, I_section, my_beam, U);
        my_beam->update_state();
    }
        void TearDown() override {
        // Reset the variables to their initial state    
        U.setZero();
        my_beam->set_global_U(U);
        my_beam->update_state();
}
};


TEST_F(PlasticBeamTests, CheckLengthCalc) 
{
    real beam_length = my_beam->get_L0();
    EXPECT_NEAR(beam_length, PLASTIC_BEAM_LENGTH, PLASTIC_BEAM_TOLERANCE);
}



TEST_F(PlasticBeamTests, RigidMoveRightCheckLocald) {
  move_right(in_nodes);
  my_beam->update_state();

  // Calculate norms and perform assertions
  vec local_d = my_beam->get_local_d();

  EXPECT_NEAR(local_d(0), 0.0, PLASTIC_BEAM_TOLERANCE);
  EXPECT_NEAR(local_d(1), 0.0, PLASTIC_BEAM_TOLERANCE);
  EXPECT_NEAR(local_d(2), 0.0, PLASTIC_BEAM_TOLERANCE);
}

TEST_F(PlasticBeamTests, RigidMoveRightCheckEps) {
  move_right(in_nodes);
  my_beam->update_state();

  // Calculate norms and perform assertions
  real eps_norm = my_beam->get_eps().lpNorm<1>();
  EXPECT_NEAR(eps_norm, 0.0, PLASTIC_BEAM_TOLERANCE);
}

TEST_F(PlasticBeamTests, RigidMoveRightCheckStress) {
  move_right(in_nodes);
  my_beam->update_state();

  // Calculate norms and perform assertions
  real stress_norm = my_beam->get_local_stresses().lpNorm<1>();
  EXPECT_NEAR(stress_norm, 0.0, PLASTIC_BEAM_TOLERANCE);
}

TEST_F(PlasticBeamTests, RigidMoveRightCheckLocalf) {
  move_right(in_nodes);
  my_beam->update_state();

  // Calculate norms and perform assertions
  real f_norm = my_beam->get_local_f().lpNorm<1>();
  EXPECT_NEAR(f_norm, 0.0, PLASTIC_BEAM_TOLERANCE);
  }

TEST_F(PlasticBeamTests, RigidMoveRightCheckResistanceForces) {
  move_right(in_nodes);
  my_beam->update_state();

  // Calculate norms and perform assertions
  real resistance_forces_norm = my_beam->get_element_resistance_forces().lpNorm<1>();
  EXPECT_NEAR(resistance_forces_norm, 0.0, PLASTIC_BEAM_TOLERANCE);
}


TEST_F(PlasticBeamTests, RigidMoveUpCheckLocald) {
  move_up(in_nodes);
  my_beam->update_state();

  // Calculate norms and perform assertions
  vec local_d = my_beam->get_local_d();

  EXPECT_NEAR(local_d(0), 0.0, PLASTIC_BEAM_TOLERANCE);
  EXPECT_NEAR(local_d(1), 0.0, PLASTIC_BEAM_TOLERANCE);
  EXPECT_NEAR(local_d(2), 0.0, PLASTIC_BEAM_TOLERANCE);
}

TEST_F(PlasticBeamTests, RigidMoveUpCheckEps) {
  move_up(in_nodes);
  my_beam->update_state();

  // Calculate norms and perform assertions
  real eps_norm = my_beam->get_eps().lpNorm<1>();
  EXPECT_NEAR(eps_norm, 0.0, PLASTIC_BEAM_TOLERANCE);
}

TEST_F(PlasticBeamTests, RigidMoveUpCheckStress) {
  move_up(in_nodes);
  my_beam->update_state();

  // Calculate norms and perform assertions
  real stress_norm = my_beam->get_local_stresses().lpNorm<1>();
  EXPECT_NEAR(stress_norm, 0.0, PLASTIC_BEAM_TOLERANCE);
}

TEST_F(PlasticBeamTests, RigidMoveUpCheckLocalf) {
  move_up(in_nodes);
  my_beam->update_state();

  // Calculate norms and perform assertions
  real f_norm = my_beam->get_local_f().lpNorm<1>();
  EXPECT_NEAR(f_norm, 0.0, PLASTIC_BEAM_TOLERANCE);
  }

TEST_F(PlasticBeamTests, RigidMoveUpCheckResistanceForces) {
  move_up(in_nodes);
  my_beam->update_state();

  // Calculate norms and perform assertions
  real resistance_forces_norm = my_beam->get_element_resistance_forces().lpNorm<1>();
  EXPECT_NEAR(resistance_forces_norm, 0.0, PLASTIC_BEAM_TOLERANCE);
}


TEST_F(PlasticBeamTests, RigidRotateCCWCheckLocald) {
  real theta = 1.0;
  real correct_extension = 0;
  std::pair<real, real> deltas = rotate_ccw(in_nodes, theta);
  real correct_thetas = theta - std::atan(2*deltas.second/(PLASTIC_BEAM_LENGTH - 2*deltas.first));
  my_beam->update_state();

  
  // Calculate norms and perform assertions
  vec local_d = my_beam->get_local_d();

  EXPECT_NEAR(local_d(0), correct_extension, PLASTIC_BEAM_TOLERANCE);
  EXPECT_NEAR(local_d(1), correct_thetas, PLASTIC_BEAM_TOLERANCE);
  EXPECT_NEAR(local_d(2), correct_thetas, PLASTIC_BEAM_TOLERANCE);
}

TEST_F(PlasticBeamTests, RigidRotateCCWCheckEps) {
  std::pair<real, real> deltas = rotate_ccw(in_nodes, 1.0);
  my_beam->update_state();

  // Calculate norms and perform assertions
  real eps_norm = my_beam->get_eps().lpNorm<1>();
  EXPECT_NEAR(eps_norm, 0.0, PLASTIC_BEAM_TOLERANCE);
}

TEST_F(PlasticBeamTests, RigidRotateCCWCheckStress) {
  std::pair<real, real> deltas = rotate_ccw(in_nodes, 1.0);
  my_beam->update_state();

  // Calculate norms and perform assertions
  real stress_norm = my_beam->get_local_stresses().lpNorm<1>();
  EXPECT_NEAR(stress_norm, 0.0, PLASTIC_BEAM_TOLERANCE);
}

TEST_F(PlasticBeamTests, RigidRotateCCWCheckLocalf) {
  std::pair<real, real> deltas = rotate_ccw(in_nodes, 1.0);
  my_beam->update_state();

  // Calculate norms and perform assertions
  real f_norm = my_beam->get_local_f().lpNorm<1>();
  EXPECT_NEAR(f_norm, 0.0, PLASTIC_BEAM_TOLERANCE);
  }

TEST_F(PlasticBeamTests, RigidRotateCCWResistanceForces) {
  std::pair<real, real> deltas = rotate_ccw(in_nodes, 1.0);
  my_beam->update_state();

  vec R = my_beam->get_element_resistance_forces();
  
  EXPECT_NEAR(R(0), 0.0, PLASTIC_BEAM_TOLERANCE); // axial
  EXPECT_NEAR(R(2), 0.0, PLASTIC_BEAM_TOLERANCE); // vertical
  EXPECT_NEAR(R(5), 0.0, PLASTIC_BEAM_TOLERANCE); // moment

  EXPECT_NEAR(R(6), 0.0, PLASTIC_BEAM_TOLERANCE); // axial
  EXPECT_NEAR(R(8), 0.0, PLASTIC_BEAM_TOLERANCE); // vertical
  EXPECT_NEAR(R(11), 0.0, PLASTIC_BEAM_TOLERANCE); // moment
}

TEST_F(PlasticBeamTests, ConstantCompressionLength) {
  real delta = PLASTIC_BEAM_LENGTH*PLASTIC_YIELD_STRENGTH_BEAM_TEST/PLASTIC_YOUNGS_MODULUS_BEAM_TEST;
  constant_compression(in_nodes, delta);
  my_beam->update_state();

  // Calculate norms and perform assertions
  real deformed_length = my_beam->get_L();
  EXPECT_NEAR(deformed_length, PLASTIC_BEAM_LENGTH - delta, PLASTIC_BEAM_TOLERANCE);
}

TEST_F(PlasticBeamTests, ConstantCompressionEps) {
  real delta = PLASTIC_BEAM_LENGTH*PLASTIC_YIELD_STRENGTH_BEAM_TEST/PLASTIC_YOUNGS_MODULUS_BEAM_TEST;
  constant_compression(in_nodes, delta);
  my_beam->update_state();

  // Calculate norms and perform assertions
  vec eps = my_beam->get_eps();
  EXPECT_NEAR(eps(0), -delta/PLASTIC_BEAM_LENGTH, PLASTIC_BEAM_TOLERANCE);
  EXPECT_NEAR(eps(1), 0.0, PLASTIC_BEAM_TOLERANCE);
}

TEST_F(PlasticBeamTests, ConstantCompressionStress) {
  real delta = PLASTIC_BEAM_LENGTH*PLASTIC_YIELD_STRENGTH_BEAM_TEST/PLASTIC_YOUNGS_MODULUS_BEAM_TEST;
  constant_compression(in_nodes, delta);
  my_beam->update_state();

  // Calculate norms and perform assertions
  vec stress = my_beam->get_local_stresses();
  EXPECT_NEAR(stress(0), -PLASTIC_YIELD_STRENGTH_BEAM_TEST*correct_area, PLASTIC_BEAM_TOLERANCE);
  EXPECT_NEAR(stress(1), 0.0, PLASTIC_BEAM_TOLERANCE);
}

TEST_F(PlasticBeamTests, ConstantCompressionLocalNodalForces) {
  real delta = PLASTIC_BEAM_LENGTH*PLASTIC_YIELD_STRENGTH_BEAM_TEST/PLASTIC_YOUNGS_MODULUS_BEAM_TEST;
  constant_compression(in_nodes, delta);
  my_beam->update_state();

  // Calculate norms and perform assertions
  vec local_f = my_beam->get_local_f();
  real local_f_norm = local_f.lpNorm<1>();
  EXPECT_NEAR(local_f(0), PLASTIC_YIELD_STRENGTH_BEAM_TEST*correct_area, PLASTIC_BEAM_TOLERANCE);
  EXPECT_NEAR(local_f(1), 0.0, PLASTIC_BEAM_TOLERANCE);
  EXPECT_NEAR(local_f(2), 0.0, PLASTIC_BEAM_TOLERANCE);
}

TEST_F(PlasticBeamTests, ConstantCompressionGlobalNodalForces) {
  real delta = PLASTIC_BEAM_LENGTH*PLASTIC_YIELD_STRENGTH_BEAM_TEST/PLASTIC_YOUNGS_MODULUS_BEAM_TEST;
  constant_compression(in_nodes, delta);
  my_beam->update_state();

  // Calculate norms and perform assertions
  vec R = my_beam->get_element_resistance_forces();
  real R_norm = R.lpNorm<1>();
  EXPECT_NEAR(R(0), -((delta/PLASTIC_BEAM_LENGTH)*PLASTIC_YOUNGS_MODULUS_BEAM_TEST*correct_area), PLASTIC_BEAM_TOLERANCE);
  EXPECT_NEAR(R(6), ((delta/PLASTIC_BEAM_LENGTH)*PLASTIC_YOUNGS_MODULUS_BEAM_TEST*correct_area), PLASTIC_BEAM_TOLERANCE);
  // the following test ensures only axial are the only forces via checking the norm.
  EXPECT_NEAR(R_norm, ((delta/PLASTIC_BEAM_LENGTH)*PLASTIC_YOUNGS_MODULUS_BEAM_TEST*correct_area)*2.0, PLASTIC_BEAM_TOLERANCE); 
}

TEST_F(PlasticBeamTests, ConstantTensionLength) {
  real delta = PLASTIC_BEAM_LENGTH*PLASTIC_YIELD_STRENGTH_BEAM_TEST/PLASTIC_YOUNGS_MODULUS_BEAM_TEST;
  constant_tension(in_nodes, delta);
  my_beam->update_state();

  // Calculate norms and perform assertions
  real deformed_length = my_beam->get_L();
  EXPECT_NEAR(deformed_length, PLASTIC_BEAM_LENGTH + delta, PLASTIC_BEAM_TOLERANCE);
}


TEST_F(PlasticBeamTests, ConstantTensionEps) {
  real delta = PLASTIC_BEAM_LENGTH*PLASTIC_YIELD_STRENGTH_BEAM_TEST/PLASTIC_YOUNGS_MODULUS_BEAM_TEST;
  constant_tension(in_nodes, delta);
  my_beam->update_state();

  // Calculate norms and perform assertions
  vec eps = my_beam->get_eps();
  EXPECT_NEAR(eps(0), (delta/PLASTIC_BEAM_LENGTH), PLASTIC_BEAM_TOLERANCE);
  EXPECT_NEAR(eps(1), 0.0, PLASTIC_BEAM_TOLERANCE);
}

TEST_F(PlasticBeamTests, ConstantTensionStress) {
  real delta = PLASTIC_BEAM_LENGTH*PLASTIC_YIELD_STRENGTH_BEAM_TEST/PLASTIC_YOUNGS_MODULUS_BEAM_TEST;
  constant_tension(in_nodes, delta);
  my_beam->update_state();

  // Calculate norms and perform assertions
  vec stress = my_beam->get_local_stresses();
  EXPECT_NEAR(stress(0), (delta/PLASTIC_BEAM_LENGTH)*PLASTIC_YOUNGS_MODULUS_BEAM_TEST*correct_area, PLASTIC_BEAM_TOLERANCE);
  EXPECT_NEAR(stress(1), 0.0, PLASTIC_BEAM_TOLERANCE);
}

TEST_F(PlasticBeamTests, ConstantTensionLocalNodalForces) {
  real delta = PLASTIC_BEAM_LENGTH*PLASTIC_YIELD_STRENGTH_BEAM_TEST/PLASTIC_YOUNGS_MODULUS_BEAM_TEST;
  constant_tension(in_nodes, delta);
  my_beam->update_state();

  // Calculate norms and perform assertions
  vec local_f = my_beam->get_local_f();
  real local_f_norm = local_f.lpNorm<1>();
  EXPECT_NEAR(local_f(0), -PLASTIC_YIELD_STRENGTH_BEAM_TEST*correct_area, PLASTIC_BEAM_TOLERANCE);
  EXPECT_NEAR(local_f(1), 0.0, PLASTIC_BEAM_TOLERANCE);
  EXPECT_NEAR(local_f(2), 0.0, PLASTIC_BEAM_TOLERANCE);
}

TEST_F(PlasticBeamTests, ConstantTensionGlobalNodalForces) {
  real delta = PLASTIC_BEAM_LENGTH*PLASTIC_YIELD_STRENGTH_BEAM_TEST/PLASTIC_YOUNGS_MODULUS_BEAM_TEST;
  constant_tension(in_nodes, delta);
  my_beam->update_state();

  // Calculate norms and perform assertions
  vec R = my_beam->get_element_resistance_forces();
  real R_norm = R.lpNorm<1>();
  EXPECT_NEAR(R(0), ((delta/PLASTIC_BEAM_LENGTH)*PLASTIC_YOUNGS_MODULUS_BEAM_TEST*correct_area), PLASTIC_BEAM_TOLERANCE);
  EXPECT_NEAR(R(6), -((delta/PLASTIC_BEAM_LENGTH)*PLASTIC_YOUNGS_MODULUS_BEAM_TEST*correct_area), PLASTIC_BEAM_TOLERANCE);
  EXPECT_NEAR(R_norm, ((delta/PLASTIC_BEAM_LENGTH)*PLASTIC_YOUNGS_MODULUS_BEAM_TEST*correct_area)*2.0, PLASTIC_BEAM_TOLERANCE);
}

//   /**
//    * @brief end rotation due to constant moment is \f$\theta = \frac{ML}{2EI}\f$. So, for a fixed end-rotation, the constant moment is \f$ M = \frac{\theta 2EI}{L}\f$.
//    * 
//    */
// TEST_F(PlasticBeamTests, ConstantRotationEps) {

//   constant_positive_bending(in_nodes);
//   my_beam->update_state();

//   // Calculate norms and perform assertions
//   vec eps = my_beam->get_eps();
//   EXPECT_NEAR(eps(0), 0.0, PLASTIC_BEAM_TOLERANCE);
//   EXPECT_NEAR(eps(1), (2/BEAM_LENGTH), PLASTIC_BEAM_TOLERANCE);
// }

// TEST_F(PlasticBeamTests, ConstantRotationStress) {
//   constant_positive_bending(in_nodes);
//   my_beam->update_state();

//   // Calculate norms and perform assertions
//   vec stress = my_beam->get_local_stresses();
//   EXPECT_NEAR(stress(0), 0.0, PLASTIC_BEAM_TOLERANCE);
//   EXPECT_NEAR(stress(1), (2*YOUNGS_MODULUS_BEAM_TEST*I/BEAM_LENGTH), PLASTIC_BEAM_TOLERANCE);
// }

// TEST_F(PlasticBeamTests, ConstantRotationLocalNodalForces) {
//   constant_positive_bending(in_nodes);
//   my_beam->update_state();

//   // Calculate norms and perform assertions
//   vec local_f = my_beam->get_local_f();
//   real local_f_norm = local_f.lpNorm<1>();
//   EXPECT_NEAR(local_f(2), -(2*YOUNGS_MODULUS_BEAM_TEST*I/BEAM_LENGTH), PLASTIC_BEAM_TOLERANCE);
//   EXPECT_NEAR(local_f(5), (2*YOUNGS_MODULUS_BEAM_TEST*I/BEAM_LENGTH), PLASTIC_BEAM_TOLERANCE);
//   EXPECT_NEAR(local_f_norm, (2*YOUNGS_MODULUS_BEAM_TEST*I/BEAM_LENGTH)*2.0, PLASTIC_BEAM_TOLERANCE);
// }

// TEST_F(PlasticBeamTests, ConstantRotationGlobalNodalForces) {
//   constant_positive_bending(in_nodes);
//   my_beam->update_state();

//   // Calculate norms and perform assertions
//   vec R = my_beam->get_element_resistance_forces();
//   real R_norm = R.lpNorm<1>();
//   EXPECT_NEAR(R(5), -(2*YOUNGS_MODULUS_BEAM_TEST*I/BEAM_LENGTH), PLASTIC_BEAM_TOLERANCE);
//   EXPECT_NEAR(R(11), (2*YOUNGS_MODULUS_BEAM_TEST*I/BEAM_LENGTH), PLASTIC_BEAM_TOLERANCE);
//   EXPECT_NEAR(R_norm, (2*YOUNGS_MODULUS_BEAM_TEST*I/BEAM_LENGTH)*2.0, PLASTIC_BEAM_TOLERANCE);
// }

#endif