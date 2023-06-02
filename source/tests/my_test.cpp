#include <iostream>
#include "gtest/gtest.h"
#include "../beam_element.hpp"
#include "../maths_defaults.hpp"
#include "../node.hpp"

#define TOLERANCE 1e-6

real get_l1_force(beam_element& my_beam, vec& d)
{
   return  (my_beam.get_K() * d).lpNorm<1>();
}

class RigidBodyMotionTest : public ::testing::Test {
    // Declare variables to be used in the fixture
public:
    std::array<node, 2> in_nodes = {node(0.0, 0.0, 0.0), node(3.0, 0.0, 0.0)};
    beam_element my_beam;
    vec d;
    void SetUp() override {
        // Create the nodes
        my_beam = beam_element(in_nodes);
        my_beam.calc_N(1.5);
        my_beam.calc_B(1.5);
        my_beam.calc_K();

        // Create the d vector
        d = make_xd_vec(6);
    }
    void TearDown() override {
        // Reset the variables to their initial state    
        d.setZero();
}
};

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