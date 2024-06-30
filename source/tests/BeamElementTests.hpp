#define TOLERANCE 1e-6
#include <iostream>
#include <vector>
#include "gtest/gtest.h"
#include "../beam_element.hpp"
#include "../BeamElementBaseClass.hpp"
#include "../BeamElementCommonInterface.hpp"
#include "../Linear2DBeamElement.hpp"
#include "../maths_defaults.hpp"
#include "../node.hpp"

real get_l1_force(std::shared_ptr<BeamElementBaseClass> my_beam, vec& d)
{
   return  (my_beam->get_k() * d).lpNorm<1>();
}

class RigidBodyMotionTest : public ::testing::Test {
    // Declare variables to be used in the fixture
public:
    std::vector<std::shared_ptr<Node>> in_nodes;
    std::shared_ptr<BeamElementBaseClass> my_beam;
    vec U;
    void SetUp() override {
        // Create the nodes
        in_nodes.push_back(std::make_shared<Node>(0.0, 0.0, 0.0));
        in_nodes.push_back(std::make_shared<Node>(3.0, 0.0, 0.0));
        
        my_beam = std::make_shared<Linear2DBeamElement>(0, in_nodes);
        // Create the d vector
        U = make_xd_vec(12);


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
    
  
    void SetUp() override {
        // Create the nodes
        in_nodes.push_back(std::make_shared<Node>(0.0, 0.0, 0.0));
        in_nodes.push_back(std::make_shared<Node>(3.0, 0.0, 0.0));
        my_beam = std::make_shared<Linear2DBeamElement>(0, in_nodes);
        
    }
    void TearDown() override {
      
}
};

class ElementStateTests : public ::testing::Test {
    // Declare variables to be used in the fixture
public:
    std::vector<std::shared_ptr<Node>> in_nodes;
    std::shared_ptr<BeamElementBaseClass> my_beam;
    vec U;
    void SetUp() override {
        // Create the nodes
        in_nodes.push_back(std::make_shared<Node>(0.0, 0.0, 0.0));
        in_nodes.push_back(std::make_shared<Node>(3.0, 0.0, 0.0));
        
        my_beam = std::make_shared<Linear2DBeamElement>(0, in_nodes);
        // Create the d vector
        U = make_xd_vec(12);


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
    EXPECT_NEAR(beam_length, 3.0, TOLERANCE);

}

TEST_F(BasicTransformationTest, CheckUnitTransformValues) {

  mat T = my_beam->get_T();
  EXPECT_NEAR(T(0,0), 1.0, TOLERANCE);
  T(0,0) = 0;
  EXPECT_NEAR(T(1,2), 1.0, TOLERANCE);
  T(1,2) = 0;
  EXPECT_NEAR(T(2,5), 1.0, TOLERANCE);
  T(2,5) = 0;
  EXPECT_NEAR(T(3,6), 1.0, TOLERANCE);
  T(3,6) = 0;
  EXPECT_NEAR(T(4,8), 1.0, TOLERANCE);
  T(4,8) = 0;
  EXPECT_NEAR(T(5,11), 1.0, TOLERANCE);
  T(5,11) = 0;

  for (int i; i < T.rows(); ++i)
  {
    for (int j; j < T.cols(); ++j)
    {
      EXPECT_NEAR(T(i,j), 0.0, TOLERANCE);
    }
  }
}


TEST_F(BasicTransformationTest, CheckOffsetUp) {
  my_beam->calc_T(0.5);

  mat T = my_beam->get_T();
  EXPECT_NEAR(T(0,5), 0.5, TOLERANCE);
  EXPECT_NEAR(T(3,11), 0.5, TOLERANCE);
}

TEST_F(BasicTransformationTest, CheckOffsetDown) {
  my_beam->calc_T(-0.5);

  mat T = my_beam->get_T();
  EXPECT_NEAR(T(0,5), -0.5, TOLERANCE);
  EXPECT_NEAR(T(3,11), -0.5, TOLERANCE);
  
}

TEST_F(BasicTransformationTest, CheckTransformedStiffnessSize) {
  mat k_g = my_beam->get_elem_global_stiffness();
  int n_cols = k_g.cols();
  int n_rows = k_g.rows();
  EXPECT_EQ(n_cols, 12);
  EXPECT_EQ(n_rows, 12);
}

TEST_F(RigidBodyMotionTest, MoveRightCheckStiffness) {
  // Modify the d vector for this test case
  U(0) = 1; // node 1 U1
  U(6) = 1; // node 2 U1

  my_beam->set_global_U(U);
  my_beam->update_state();

  // Calculate norms and perform assertions
  real eps_norm = my_beam->get_eps().lpNorm<1>();
  real f_norm = my_beam->get_local_f().lpNorm<1>();
  real stress_norm = my_beam->get_local_stresses().lpNorm<1>();
  real resistance_forces_norm = my_beam->get_element_resistance_forces().lpNorm<1>();
  EXPECT_NEAR(eps_norm, 0.0, TOLERANCE);
  EXPECT_NEAR(f_norm, 0.0, TOLERANCE);
  EXPECT_NEAR(stress_norm, 0.0, TOLERANCE);
  EXPECT_NEAR(resistance_forces_norm, 0.0, TOLERANCE);
}

TEST_F(RigidBodyMotionTest, MoveUpCheckStiffness) {
  // Modify the d vector for this test case
  U(2) = 1; // node 1 U2
  U(8) = 1; // node 2 U2
  my_beam->set_global_U(U);
  my_beam->update_state();

  // Calculate norms and perform assertions
  real eps_norm = my_beam->get_eps().lpNorm<1>();
  real f_norm = my_beam->get_local_f().lpNorm<1>();
  real stress_norm = my_beam->get_local_stresses().lpNorm<1>();
  real resistance_forces_norm = my_beam->get_element_resistance_forces().lpNorm<1>();
  EXPECT_NEAR(eps_norm, 0.0, TOLERANCE);
  EXPECT_NEAR(f_norm, 0.0, TOLERANCE);
  EXPECT_NEAR(stress_norm, 0.0, TOLERANCE);
  EXPECT_NEAR(resistance_forces_norm, 0.0, TOLERANCE);
  
}

TEST_F(RigidBodyMotionTest, RotateCCWCheckStiffness) {
  // Modify the d vector for this test case
  U(2) = -1; // node 1 U2
  U(5) = 2.0/3; // node 1 U33
  U(8) = 1; // node 2 U2
  U(11) = 2.0/3; // node 2 U33

  my_beam->set_global_U(U);
  my_beam->update_state();

  // Calculate norms and perform assertions
  real eps_norm = my_beam->get_eps().lpNorm<1>();
  real f_norm = my_beam->get_local_f().lpNorm<1>();
  real stress_norm = my_beam->get_local_stresses().lpNorm<1>();
  real resistance_forces_norm = my_beam->get_element_resistance_forces().lpNorm<1>();
  EXPECT_NEAR(eps_norm, 0.0, TOLERANCE);
  EXPECT_NEAR(f_norm, 0.0, TOLERANCE);
  EXPECT_NEAR(stress_norm, 0.0, TOLERANCE);
  EXPECT_NEAR(resistance_forces_norm, 0.0, TOLERANCE);
  
}
