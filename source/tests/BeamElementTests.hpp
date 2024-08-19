/**
 * @file BeamElementTests.hpp
 * @brief tests the linear beam-column element. Careful as some of these tests not suitable for the nonlinear beam-column.
 * 
 */
#ifndef BEAM_ELEMENT_TESTS_HPP
#define BEAM_ELEMENT_TESTS_HPP

#define BEAM_TOLERANCE 1e-6
#include <iostream>
#include <vector>
#include "gtest/gtest.h"
#include "../beam_element.hpp"
#include "../BeamElementBaseClass.hpp"
#include "../BeamElementCommonInterface.hpp"
#include "../Linear2DBeamElement.hpp"
#include "../Nonlinear2DBeamElement.hpp"
#include "../maths_defaults.hpp"
#include "../node.hpp"

#define BEAM_LENGTH 3.0
#define YOUNGS_MODULUS_BEAM_TEST 2.06e11
#define A 0.0125
#define I 0.0004570000

void common_beam_setup(std::vector<std::shared_ptr<Node>>& in_nodes, std::shared_ptr<BeamElementBaseClass>& my_beam, vec& U) 
{
    // Create the nodes
    in_nodes.push_back(std::make_shared<Node>(0.0, 0.0, 0.0));
    in_nodes.push_back(std::make_shared<Node>(BEAM_LENGTH, 0.0, 0.0));
    BasicSection sect(2.06e11, 0.0125, 0.0004570000);
    my_beam = std::make_shared<Linear2DBeamElement>(0, in_nodes, sect);
    // Create the d vector
    U = make_xd_vec(12);
}
void move_right(std::vector<std::shared_ptr<Node>>& in_nodes) 
{
    in_nodes[0]->set_nodal_displacement(0, 1.0); // node 1 U1
    in_nodes[1]->set_nodal_displacement(0, 1.0); // node 2 U1
}
void move_up(std::vector<std::shared_ptr<Node>>& in_nodes) 
{
    in_nodes[0]->set_nodal_displacement(2, 1.0); // node 1 U2
    in_nodes[1]->set_nodal_displacement(2, 1.0); // node 2 U2
}
void rotate_ccw(std::vector<std::shared_ptr<Node>>& in_nodes) 
{
    in_nodes[0]->set_nodal_displacement(2, -1.0); // node 1 U2
    in_nodes[0]->set_nodal_displacement(5, 2.0/BEAM_LENGTH); // node 1 U33
    in_nodes[1]->set_nodal_displacement(2, 1.0); // node 2 U2
    in_nodes[1]->set_nodal_displacement(5, 2.0/BEAM_LENGTH); // node 2 U33
}

void constant_compression(std::vector<std::shared_ptr<Node>>& in_nodes) 
{
    in_nodes[0]->set_nodal_displacement(0, 0.5); // node 1 U1
    in_nodes[1]->set_nodal_displacement(0, -0.5); // node 2 U1
}
void constant_tension(std::vector<std::shared_ptr<Node>>& in_nodes) 
{
    in_nodes[0]->set_nodal_displacement(0, -0.5); // node 1 U1
    in_nodes[1]->set_nodal_displacement(0, 0.5); // node 2 U1
}
void constant_positive_bending(std::vector<std::shared_ptr<Node>>& in_nodes) 
{
    in_nodes[0]->set_nodal_displacement(5, -1); // node 2 U33
    in_nodes[1]->set_nodal_displacement(5, 1); // node 2 U33
}

class RigidBodyMotionTest : public ::testing::Test {
    // Declare variables to be used in the fixture
public:
    std::vector<std::shared_ptr<Node>> in_nodes;
    std::shared_ptr<BeamElementBaseClass> my_beam;
    vec U;
    void SetUp() override {
      common_beam_setup(in_nodes, my_beam, U);
    }
    void TearDown() override {
        // Reset the variables to their initial state    
        U.setZero();
        my_beam->set_global_U(U);
        my_beam->update_state();
}
};
class BasicTransformationTest : public ::testing::Test {
    // Declare variables to be used in the fixture
  public:
    std::vector<std::shared_ptr<Node>> in_nodes;
    std::shared_ptr<BeamElementBaseClass> my_beam;
    vec U; // not needed for this test but I don't want to create another \ref common_beam_setup function.
    
  
    void SetUp() override {
        // Create the nodes
        common_beam_setup(in_nodes, my_beam, U);  
    }
    void TearDown() override {
      
}
};

class ConstantStrainStateTest : public ::testing::Test {
    // Declare variables to be used in the fixture
  public:
    std::vector<std::shared_ptr<Node>> in_nodes;
    std::shared_ptr<BeamElementBaseClass> my_beam;
    vec U;
    void SetUp() override {
        common_beam_setup(in_nodes, my_beam, U);
    }
    void TearDown() override {
      // Reset the variables to their initial state    
      U.setZero();
      my_beam->set_global_U(U);
      my_beam->update_state();
    }
};

class ElementMappingTest : public ::testing::Test {
  public:
    std::vector<std::shared_ptr<Node>> in_nodes;
    std::shared_ptr<BeamElementBaseClass> my_beam;
    vec U;
    void SetUp() override {
      common_beam_setup(in_nodes, my_beam, U);
    }

    void TearDown() override {
      // Reset the variables to their initial state    
      U.setZero();
      my_beam->set_global_U(U);
      my_beam->update_state();
    }
};

TEST_F(BasicTransformationTest, CheckLengthCalc) 
{
    real beam_length = my_beam->get_L();
    EXPECT_NEAR(beam_length, 3.0, BEAM_TOLERANCE);

}

TEST_F(BasicTransformationTest, CheckUnitTransformValues) {

  mat T = my_beam->get_T();
  EXPECT_NEAR(T(0,0), 1.0, BEAM_TOLERANCE);
  T(0,0) = 0;
  EXPECT_NEAR(T(1,2), 1.0, BEAM_TOLERANCE);
  T(1,2) = 0;
  EXPECT_NEAR(T(2,5), 1.0, BEAM_TOLERANCE);
  T(2,5) = 0;
  EXPECT_NEAR(T(3,6), 1.0, BEAM_TOLERANCE);
  T(3,6) = 0;
  EXPECT_NEAR(T(4,8), 1.0, BEAM_TOLERANCE);
  T(4,8) = 0;
  EXPECT_NEAR(T(5,11), 1.0, BEAM_TOLERANCE);
  T(5,11) = 0;

  for (int i; i < T.rows(); ++i)
  {
    for (int j; j < T.cols(); ++j)
    {
      EXPECT_NEAR(T(i,j), 0.0, BEAM_TOLERANCE);
    }
  }
}


TEST_F(BasicTransformationTest, CheckOffsetUp) {
  my_beam->calc_T(0.5);

  mat T = my_beam->get_T();
  EXPECT_NEAR(T(0,5), 0.5, BEAM_TOLERANCE);
  EXPECT_NEAR(T(3,11), 0.5, BEAM_TOLERANCE);
}

TEST_F(BasicTransformationTest, CheckOffsetDown) {
  my_beam->calc_T(-0.5);

  mat T = my_beam->get_T();
  EXPECT_NEAR(T(0,5), -0.5, BEAM_TOLERANCE);
  EXPECT_NEAR(T(3,11), -0.5, BEAM_TOLERANCE);
  
}

TEST_F(BasicTransformationTest, CheckTransformedStiffnessSize) {
  mat k_g = my_beam->get_elem_global_stiffness();
  int n_cols = k_g.cols();
  int n_rows = k_g.rows();
  EXPECT_EQ(n_cols, 12);
  EXPECT_EQ(n_rows, 12);
}


TEST_F(RigidBodyMotionTest, MoveRightCheckLocald) {
  move_right(in_nodes);
  my_beam->update_state();

  // Calculate norms and perform assertions
  vec local_d = my_beam->get_local_d();

  EXPECT_NEAR(local_d(0), 1.0, BEAM_TOLERANCE);
  EXPECT_NEAR(local_d(1), 0.0, BEAM_TOLERANCE);
  EXPECT_NEAR(local_d(2), 0.0, BEAM_TOLERANCE);
  EXPECT_NEAR(local_d(3), 1.0, BEAM_TOLERANCE);
  EXPECT_NEAR(local_d(4), 0.0, BEAM_TOLERANCE);
  EXPECT_NEAR(local_d(5), 0.0, BEAM_TOLERANCE);
}

TEST_F(RigidBodyMotionTest, MoveRightCheckEps) {
  move_right(in_nodes);
  my_beam->update_state();

  // Calculate norms and perform assertions
  real eps_norm = my_beam->get_eps().lpNorm<1>();
  EXPECT_NEAR(eps_norm, 0.0, BEAM_TOLERANCE);
}

TEST_F(RigidBodyMotionTest, MoveRightCheckStress) {
  move_right(in_nodes);
  my_beam->update_state();

  // Calculate norms and perform assertions
  real stress_norm = my_beam->get_local_stresses().lpNorm<1>();
  EXPECT_NEAR(stress_norm, 0.0, BEAM_TOLERANCE);
}

TEST_F(RigidBodyMotionTest, MoveRightCheckLocalf) {
  move_right(in_nodes);
  my_beam->update_state();

  // Calculate norms and perform assertions
  real f_norm = my_beam->get_local_f().lpNorm<1>();
  EXPECT_NEAR(f_norm, 0.0, BEAM_TOLERANCE);
  }

TEST_F(RigidBodyMotionTest, MoveRightCheckResistanceForces) {
  move_right(in_nodes);
  my_beam->update_state();

  // Calculate norms and perform assertions
  real resistance_forces_norm = my_beam->get_element_resistance_forces().lpNorm<1>();
  EXPECT_NEAR(resistance_forces_norm, 0.0, BEAM_TOLERANCE);
}


TEST_F(RigidBodyMotionTest, MoveUpCheckLocald) {
  move_up(in_nodes);
  my_beam->update_state();

  // Calculate norms and perform assertions
  vec local_d = my_beam->get_local_d();

  EXPECT_NEAR(local_d(0), 0.0, BEAM_TOLERANCE);
  EXPECT_NEAR(local_d(1), 1.0, BEAM_TOLERANCE);
  EXPECT_NEAR(local_d(2), 0.0, BEAM_TOLERANCE);
  EXPECT_NEAR(local_d(3), 0.0, BEAM_TOLERANCE);
  EXPECT_NEAR(local_d(4), 1.0, BEAM_TOLERANCE);
  EXPECT_NEAR(local_d(5), 0.0, BEAM_TOLERANCE);
}

TEST_F(RigidBodyMotionTest, MoveUpCheckEps) {
  move_up(in_nodes);
  my_beam->update_state();

  // Calculate norms and perform assertions
  real eps_norm = my_beam->get_eps().lpNorm<1>();
  EXPECT_NEAR(eps_norm, 0.0, BEAM_TOLERANCE);
}

TEST_F(RigidBodyMotionTest, MoveUpCheckStress) {
  move_up(in_nodes);
  my_beam->update_state();

  // Calculate norms and perform assertions
  real stress_norm = my_beam->get_local_stresses().lpNorm<1>();
  EXPECT_NEAR(stress_norm, 0.0, BEAM_TOLERANCE);
}

TEST_F(RigidBodyMotionTest, MoveUpCheckLocalf) {
  move_up(in_nodes);
  my_beam->update_state();

  // Calculate norms and perform assertions
  real f_norm = my_beam->get_local_f().lpNorm<1>();
  EXPECT_NEAR(f_norm, 0.0, BEAM_TOLERANCE);
  }

TEST_F(RigidBodyMotionTest, MoveUpCheckResistanceForces) {
  move_up(in_nodes);
  my_beam->update_state();

  // Calculate norms and perform assertions
  real resistance_forces_norm = my_beam->get_element_resistance_forces().lpNorm<1>();
  EXPECT_NEAR(resistance_forces_norm, 0.0, BEAM_TOLERANCE);
}


TEST_F(RigidBodyMotionTest, RotateCCWCheckLocald) {
  rotate_ccw(in_nodes);
  my_beam->update_state();

  // Calculate norms and perform assertions
  vec local_d = my_beam->get_local_d();

  EXPECT_NEAR(local_d(0), 0.0, BEAM_TOLERANCE);
  EXPECT_NEAR(local_d(1), -1.0, BEAM_TOLERANCE);
  EXPECT_NEAR(local_d(2), (2.0/BEAM_LENGTH), BEAM_TOLERANCE);
  EXPECT_NEAR(local_d(3), 0.0, BEAM_TOLERANCE);
  EXPECT_NEAR(local_d(4), 1.0, BEAM_TOLERANCE);
  EXPECT_NEAR(local_d(5), (2.0/BEAM_LENGTH), BEAM_TOLERANCE);
}

TEST_F(RigidBodyMotionTest, RotateCCWCheckEps) {
  rotate_ccw(in_nodes);
  my_beam->update_state();

  // Calculate norms and perform assertions
  real eps_norm = my_beam->get_eps().lpNorm<1>();
  EXPECT_NEAR(eps_norm, 0.0, BEAM_TOLERANCE);
}

TEST_F(RigidBodyMotionTest, RotateCCWCheckStress) {
  rotate_ccw(in_nodes);
  my_beam->update_state();

  // Calculate norms and perform assertions
  real stress_norm = my_beam->get_local_stresses().lpNorm<1>();
  EXPECT_NEAR(stress_norm, 0.0, BEAM_TOLERANCE);
}

TEST_F(RigidBodyMotionTest, RotateCCWCheckLocalf) {
  rotate_ccw(in_nodes);
  my_beam->update_state();

  // Calculate norms and perform assertions
  real f_norm = my_beam->get_local_f().lpNorm<1>();
  EXPECT_NEAR(f_norm, 0.0, BEAM_TOLERANCE);
  }

TEST_F(RigidBodyMotionTest, RotateCCWResistanceForces) {
  rotate_ccw(in_nodes);
  my_beam->update_state();

  // Calculate norms and perform assertions
  real resistance_forces_norm = my_beam->get_element_resistance_forces().lpNorm<1>();
  EXPECT_NEAR(resistance_forces_norm, 0.0, BEAM_TOLERANCE);
}


TEST_F(ConstantStrainStateTest, ConstantCompressionEps) {
  constant_compression(in_nodes);
  my_beam->update_state();

  // Calculate norms and perform assertions
  vec eps = my_beam->get_eps();
  EXPECT_NEAR(eps(0), (-1.0/BEAM_LENGTH), BEAM_TOLERANCE);
  EXPECT_NEAR(eps(1), 0.0, BEAM_TOLERANCE);
}

TEST_F(ConstantStrainStateTest, ConstantCompressionStress) {
  constant_compression(in_nodes);
  my_beam->update_state();

  // Calculate norms and perform assertions
  vec stress = my_beam->get_local_stresses();
  EXPECT_NEAR(stress(0), (-1.0/BEAM_LENGTH)*YOUNGS_MODULUS_BEAM_TEST*A, BEAM_TOLERANCE);
  EXPECT_NEAR(stress(1), 0.0, BEAM_TOLERANCE);
}

TEST_F(ConstantStrainStateTest, ConstantCompressionLocalNodalForces) {
  constant_compression(in_nodes);
  my_beam->update_state();

  // Calculate norms and perform assertions
  vec local_f = my_beam->get_local_f();
  real local_f_norm = local_f.lpNorm<1>();
  EXPECT_NEAR(local_f(0), ((1.0/BEAM_LENGTH)*YOUNGS_MODULUS_BEAM_TEST*A), BEAM_TOLERANCE);
  EXPECT_NEAR(local_f(3), -((1.0/BEAM_LENGTH)*YOUNGS_MODULUS_BEAM_TEST*A), BEAM_TOLERANCE);
  EXPECT_NEAR(local_f_norm, ((1.0/BEAM_LENGTH)*YOUNGS_MODULUS_BEAM_TEST*A)*2.0, BEAM_TOLERANCE);
}

TEST_F(ConstantStrainStateTest, ConstantCompressionGlobalNodalForces) {
  constant_compression(in_nodes);
  my_beam->update_state();

  // Calculate norms and perform assertions
  vec R = my_beam->get_element_resistance_forces();
  real R_norm = R.lpNorm<1>();
  EXPECT_NEAR(R(0), ((1.0/BEAM_LENGTH)*YOUNGS_MODULUS_BEAM_TEST*A), BEAM_TOLERANCE);
  EXPECT_NEAR(R(6), -((1.0/BEAM_LENGTH)*YOUNGS_MODULUS_BEAM_TEST*A), BEAM_TOLERANCE);
  EXPECT_NEAR(R_norm, ((1.0/BEAM_LENGTH)*YOUNGS_MODULUS_BEAM_TEST*A)*2.0, BEAM_TOLERANCE);
}

TEST_F(ConstantStrainStateTest, ConstantTensionEps) {
  constant_tension(in_nodes);
  my_beam->update_state();

  // Calculate norms and perform assertions
  vec eps = my_beam->get_eps();
  EXPECT_NEAR(eps(0), (1.0/BEAM_LENGTH), BEAM_TOLERANCE);
  EXPECT_NEAR(eps(1), 0.0, BEAM_TOLERANCE);
}

TEST_F(ConstantStrainStateTest, ConstantTensionStress) {
  constant_tension(in_nodes);
  my_beam->update_state();

  // Calculate norms and perform assertions
  vec stress = my_beam->get_local_stresses();
  EXPECT_NEAR(stress(0), (1.0/BEAM_LENGTH)*YOUNGS_MODULUS_BEAM_TEST*A, BEAM_TOLERANCE);
  EXPECT_NEAR(stress(1), 0.0, BEAM_TOLERANCE);
}

TEST_F(ConstantStrainStateTest, ConstantTensionLocalNodalForces) {
  constant_tension(in_nodes);
  my_beam->update_state();

  // Calculate norms and perform assertions
  vec local_f = my_beam->get_local_f();
  real local_f_norm = local_f.lpNorm<1>();
  EXPECT_NEAR(local_f(0), -((1.0/BEAM_LENGTH)*YOUNGS_MODULUS_BEAM_TEST*A), BEAM_TOLERANCE);
  EXPECT_NEAR(local_f(3), ((1.0/BEAM_LENGTH)*YOUNGS_MODULUS_BEAM_TEST*A), BEAM_TOLERANCE);
  EXPECT_NEAR(local_f_norm, ((1.0/BEAM_LENGTH)*YOUNGS_MODULUS_BEAM_TEST*A)*2.0, BEAM_TOLERANCE);
}

TEST_F(ConstantStrainStateTest, ConstantTensionGlobalNodalForces) {
  constant_tension(in_nodes);
  my_beam->update_state();

  // Calculate norms and perform assertions
  vec R = my_beam->get_element_resistance_forces();
  real R_norm = R.lpNorm<1>();
  EXPECT_NEAR(R(0), -((1.0/BEAM_LENGTH)*YOUNGS_MODULUS_BEAM_TEST*A), BEAM_TOLERANCE);
  EXPECT_NEAR(R(6), ((1.0/BEAM_LENGTH)*YOUNGS_MODULUS_BEAM_TEST*A), BEAM_TOLERANCE);
  EXPECT_NEAR(R_norm, ((1.0/BEAM_LENGTH)*YOUNGS_MODULUS_BEAM_TEST*A)*2.0, BEAM_TOLERANCE);
}

  /**
   * @brief end rotation due to constant moment is \f$\theta = \frac{ML}{2EI}\f$. So, for a fixed end-rotation, the constant moment is \f$ M = \frac{\theta 2EI}{L}\f$.
   * 
   */
TEST_F(ConstantStrainStateTest, ConstantRotationEps) {

  constant_positive_bending(in_nodes);
  my_beam->update_state();

  // Calculate norms and perform assertions
  vec eps = my_beam->get_eps();
  EXPECT_NEAR(eps(0), 0.0, BEAM_TOLERANCE);
  EXPECT_NEAR(eps(1), (2/BEAM_LENGTH), BEAM_TOLERANCE);
}

TEST_F(ConstantStrainStateTest, ConstantRotationStress) {
  constant_positive_bending(in_nodes);
  my_beam->update_state();

  // Calculate norms and perform assertions
  vec stress = my_beam->get_local_stresses();
  EXPECT_NEAR(stress(0), 0.0, BEAM_TOLERANCE);
  EXPECT_NEAR(stress(1), (2*YOUNGS_MODULUS_BEAM_TEST*I/BEAM_LENGTH), BEAM_TOLERANCE);
}

TEST_F(ConstantStrainStateTest, ConstantRotationLocalNodalForces) {
  constant_positive_bending(in_nodes);
  my_beam->update_state();

  // Calculate norms and perform assertions
  vec local_f = my_beam->get_local_f();
  real local_f_norm = local_f.lpNorm<1>();
  EXPECT_NEAR(local_f(2), -(2*YOUNGS_MODULUS_BEAM_TEST*I/BEAM_LENGTH), BEAM_TOLERANCE);
  EXPECT_NEAR(local_f(5), (2*YOUNGS_MODULUS_BEAM_TEST*I/BEAM_LENGTH), BEAM_TOLERANCE);
  EXPECT_NEAR(local_f_norm, (2*YOUNGS_MODULUS_BEAM_TEST*I/BEAM_LENGTH)*2.0, BEAM_TOLERANCE);
}

TEST_F(ConstantStrainStateTest, ConstantRotationGlobalNodalForces) {
  constant_positive_bending(in_nodes);
  my_beam->update_state();

  // Calculate norms and perform assertions
  vec R = my_beam->get_element_resistance_forces();
  real R_norm = R.lpNorm<1>();
  EXPECT_NEAR(R(5), -(2*YOUNGS_MODULUS_BEAM_TEST*I/BEAM_LENGTH), BEAM_TOLERANCE);
  EXPECT_NEAR(R(11), (2*YOUNGS_MODULUS_BEAM_TEST*I/BEAM_LENGTH), BEAM_TOLERANCE);
  EXPECT_NEAR(R_norm, (2*YOUNGS_MODULUS_BEAM_TEST*I/BEAM_LENGTH)*2.0, BEAM_TOLERANCE);
}

#endif