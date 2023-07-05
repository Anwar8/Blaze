#include <iostream>
#include "gtest/gtest.h"
#include "../beam_element.hpp"
#include "../maths_defaults.hpp"
#include "../node.hpp"

#define TOLERANCE 1e-6

real get_l1_force(Basic2DBeamElement& my_beam, vec& d)
{
   return  (my_beam.get_k() * d).lpNorm<1>();
}

class RigidBodyMotionTest : public ::testing::Test {
    // Declare variables to be used in the fixture
public:
    
    std::shared_ptr<Node> in_nodes_1 = std::make_shared<Node>(0.0, 0.0, 0.0);
    std::shared_ptr<Node> in_nodes_2 = std::make_shared<Node>(3.0, 0.0, 0.0);
    Basic2DBeamElement my_beam;
    vec d;
    void SetUp() override {
        // Create the nodes
        my_beam = Basic2DBeamElement(in_nodes_1, in_nodes_2);
        my_beam.calc_N(1.5);
        my_beam.calc_B(1.5);
        my_beam.calc_k();

        // Create the d vector
        d = make_xd_vec(6);
    }
    void TearDown() override {
        // Reset the variables to their initial state    
        d.setZero();
}
};
class BasicTransformationTest : public ::testing::Test {
    // Declare variables to be used in the fixture
public:
    
    std::shared_ptr<Node> in_nodes_1 = std::make_shared<Node>(0.0, 0.0, 0.0);
    std::shared_ptr<Node> in_nodes_2 = std::make_shared<Node>(3.0, 0.0, 0.0);
    Basic2DBeamElement my_beam;
    void SetUp() override {
        // Create the nodes
        my_beam = Basic2DBeamElement(in_nodes_1, in_nodes_2);
        my_beam.calc_k();
        my_beam.calc_T();
    }
    void TearDown() override {
      
}
};

TEST_F(BasicTransformationTest, CheckUnitTransformValues) {

  mat T = my_beam.get_T();
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
  my_beam.calc_T(0.5);

  mat T = my_beam.get_T();
  EXPECT_NEAR(T(0,5), 0.5, TOLERANCE);
  EXPECT_NEAR(T(1,5), 0.0, TOLERANCE);
  EXPECT_NEAR(T(3,11), 0.5, TOLERANCE);
  EXPECT_NEAR(T(4,11), 0.0, TOLERANCE);
}

TEST_F(BasicTransformationTest, CheckOffsetDown) {
  my_beam.calc_T(-0.5);

  mat T = my_beam.get_T();
  EXPECT_NEAR(T(0,5), -0.5, TOLERANCE);
  EXPECT_NEAR(T(1,5), 0.0, TOLERANCE);
  EXPECT_NEAR(T(3,11), -0.5, TOLERANCE);
  EXPECT_NEAR(T(4,11), 0.0, TOLERANCE);
}

TEST_F(BasicTransformationTest, CheckTransformedStiffnessSize) {

  mat T = my_beam.get_T();
  mat k = my_beam.get_k();
  mat k_g = T.transpose()*k*T;
  int n_cols = k_g.cols();
  int n_rows = k_g.rows();
  EXPECT_EQ(n_cols, 12);
  EXPECT_EQ(n_rows, 12);
}

TEST_F(RigidBodyMotionTest, MoveRightCheckStiffness) {
  // Modify the d vector for this test case
  d(0) = 1;
  d(3) = 1;
  my_beam.set_d(d);

  // Calculate eps and perform assertions
  
  real f = get_l1_force(my_beam, d);
  EXPECT_NEAR(f, 0.0, TOLERANCE);
  
}

TEST_F(RigidBodyMotionTest, MoveUpCheckStiffness) {
  // Modify the d vector for this test case
  d(1) = 1;
  d(4) = 1;
  my_beam.set_d(d);

  // Calculate eps and perform assertions
  
  real f = get_l1_force(my_beam, d);
  EXPECT_NEAR(f, 0.0, TOLERANCE);
  
}

TEST_F(RigidBodyMotionTest, RotateCCWCheckStiffness) {
  // Modify the d vector for this test case
  d(1) = -1; 
  d(2) = 2.0/3;
  d(4) = 1;
  d(5) = 2.0/3;
  my_beam.set_d(d);

  // Calculate eps and perform assertions
  
  real f = get_l1_force(my_beam, d);
  EXPECT_NEAR(f, 0.0, TOLERANCE);
  
}

TEST_F(RigidBodyMotionTest, MoveRight) {
  // Modify the d vector for this test case
  d(0) = 1;
  d(3) = 1;
  my_beam.set_d(d);

  // Calculate eps and perform assertions
  my_beam.calc_eps();
  vec eps = my_beam.get_eps();
  EXPECT_NEAR(eps(0), 0.0, TOLERANCE);
  EXPECT_NEAR(eps(1), 0.0, TOLERANCE);
}

TEST_F(RigidBodyMotionTest, MoveUp) {
  // Modify the d vector for this test case
  d(1) = 1;
  d(4) = 1;
  my_beam.set_d(d);

  // Calculate eps and perform assertions
  my_beam.calc_eps();
  vec eps = my_beam.get_eps();
  EXPECT_NEAR(eps(0), 0.0, TOLERANCE);
  EXPECT_NEAR(eps(1), 0.0, TOLERANCE);
}

TEST_F(RigidBodyMotionTest, RotateCCW) {
  // Modify the d vector for this test case
  d(1) = -1; 
  d(2) = 2.0/3;
  d(4) = 1;
  d(5) = 2.0/3;
  my_beam.set_d(d);

  // Calculate eps and perform assertions
  my_beam.calc_eps();
  vec eps = my_beam.get_eps();
  EXPECT_NEAR(eps(0), 0.0, TOLERANCE);
  EXPECT_NEAR(eps(1), 0.0, TOLERANCE);
  
}


int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}