/**
 * @file PlasticBeamElementTests.hpp
 * @brief tests the nonlinear plastic beam-column element.
 * 
 */
#ifndef PLASTIC_BEAM_ELEMENT_TESTS_HPP
#define PLASTIC_BEAM_ELEMENT_TESTS_HPP

#include "TestHelpers.hpp"

class PlasticBeamTests : public ::testing::Test {
  public:
    
    
    std::vector<std::shared_ptr<Node>> in_nodes;
    std::shared_ptr<Nonlinear2DPlasticBeamElement> my_beam;
    vec U;
    CommonSectionDefinitions common;
    void SetUp() override {
        common.initialise_section();
        common_plastic_beam_setup(PLASTIC_BEAM_LENGTH, in_nodes, common.I_section, my_beam, U);
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
    EXPECT_NEAR(beam_length, PLASTIC_BEAM_LENGTH, BASIC_TOLERANCE);
}


TEST_F(PlasticBeamTests, RigidMoveRightCheckLocald) {
  move_laterally(in_nodes);
  my_beam->update_state();

  // Calculate norms and perform assertions
  vec local_d = my_beam->get_local_d();

  EXPECT_NEAR(local_d(0), 0.0, BASIC_TOLERANCE);
  EXPECT_NEAR(local_d(1), 0.0, BASIC_TOLERANCE);
  EXPECT_NEAR(local_d(2), 0.0, BASIC_TOLERANCE);
}

TEST_F(PlasticBeamTests, RigidMoveRightCheckEps) {
  move_laterally(in_nodes);
  my_beam->update_state();

  // Calculate norms and perform assertions
  real eps_norm = my_beam->get_eps().lpNorm<1>();
  EXPECT_NEAR(eps_norm, 0.0, BASIC_TOLERANCE);
}

TEST_F(PlasticBeamTests, RigidMoveRightCheckStress) {
  move_laterally(in_nodes);
  my_beam->update_state();

  // Calculate norms and perform assertions
  real stress_norm = my_beam->get_local_stresses().lpNorm<1>();
  EXPECT_NEAR(stress_norm, 0.0, BASIC_TOLERANCE);
}

TEST_F(PlasticBeamTests, RigidMoveRightCheckLocalf) {
  move_laterally(in_nodes);
  my_beam->update_state();

  // Calculate norms and perform assertions
  real f_norm = my_beam->get_local_f().lpNorm<1>();
  EXPECT_NEAR(f_norm, 0.0, BASIC_TOLERANCE);
  }

TEST_F(PlasticBeamTests, RigidMoveRightCheckResistanceForces) {
  move_laterally(in_nodes);
  my_beam->update_state();

  // Calculate norms and perform assertions
  real resistance_forces_norm = my_beam->get_element_resistance_forces().lpNorm<1>();
  EXPECT_NEAR(resistance_forces_norm, 0.0, BASIC_TOLERANCE);
}


TEST_F(PlasticBeamTests, RigidMoveUpCheckLocald) {
  move_vertically(in_nodes);
  my_beam->update_state();

  // Calculate norms and perform assertions
  vec local_d = my_beam->get_local_d();

  EXPECT_NEAR(local_d(0), 0.0, BASIC_TOLERANCE);
  EXPECT_NEAR(local_d(1), 0.0, BASIC_TOLERANCE);
  EXPECT_NEAR(local_d(2), 0.0, BASIC_TOLERANCE);
}

TEST_F(PlasticBeamTests, RigidMoveUpCheckEps) {
  move_vertically(in_nodes);
  my_beam->update_state();

  // Calculate norms and perform assertions
  real eps_norm = my_beam->get_eps().lpNorm<1>();
  EXPECT_NEAR(eps_norm, 0.0, BASIC_TOLERANCE);
}

TEST_F(PlasticBeamTests, RigidMoveUpCheckStress) {
  move_vertically(in_nodes);
  my_beam->update_state();

  // Calculate norms and perform assertions
  real stress_norm = my_beam->get_local_stresses().lpNorm<1>();
  EXPECT_NEAR(stress_norm, 0.0, BASIC_TOLERANCE);
}

TEST_F(PlasticBeamTests, RigidMoveUpCheckLocalf) {
  move_vertically(in_nodes);
  my_beam->update_state();

  // Calculate norms and perform assertions
  real f_norm = my_beam->get_local_f().lpNorm<1>();
  EXPECT_NEAR(f_norm, 0.0, BASIC_TOLERANCE);
  }

TEST_F(PlasticBeamTests, RigidMoveUpCheckResistanceForces) {
  move_vertically(in_nodes);
  my_beam->update_state();

  // Calculate norms and perform assertions
  real resistance_forces_norm = my_beam->get_element_resistance_forces().lpNorm<1>();
  EXPECT_NEAR(resistance_forces_norm, 0.0, BASIC_TOLERANCE);
}


TEST_F(PlasticBeamTests, RigidRotateCCWCheckLocald) {
  real theta = 1.0;
  real correct_extension = 0;
  std::pair<real, real> deltas = rotate_ccw(in_nodes, theta);
  real correct_thetas = theta - std::atan(2*deltas.second/(PLASTIC_BEAM_LENGTH - 2*deltas.first));
  my_beam->update_state();

  
  // Calculate norms and perform assertions
  vec local_d = my_beam->get_local_d();

  EXPECT_NEAR(local_d(0), correct_extension, BASIC_TOLERANCE);
  EXPECT_NEAR(local_d(1), correct_thetas, BASIC_TOLERANCE);
  EXPECT_NEAR(local_d(2), correct_thetas, BASIC_TOLERANCE);
}

TEST_F(PlasticBeamTests, RigidRotateCCWCheckEps) {
  std::pair<real, real> deltas = rotate_ccw(in_nodes, 1.0);
  my_beam->update_state();

  // Calculate norms and perform assertions
  real eps_norm = my_beam->get_eps().lpNorm<1>();
  EXPECT_NEAR(eps_norm, 0.0, BASIC_TOLERANCE);
}

TEST_F(PlasticBeamTests, RigidRotateCCWCheckStress) {
  std::pair<real, real> deltas = rotate_ccw(in_nodes, 1.0);
  my_beam->update_state();

  // Calculate norms and perform assertions
  real stress_norm = my_beam->get_local_stresses().lpNorm<1>();
  EXPECT_NEAR(stress_norm, 0.0, BASIC_TOLERANCE);
}

TEST_F(PlasticBeamTests, RigidRotateCCWCheckLocalf) {
  std::pair<real, real> deltas = rotate_ccw(in_nodes, 1.0);
  my_beam->update_state();

  // Calculate norms and perform assertions
  real f_norm = my_beam->get_local_f().lpNorm<1>();
  EXPECT_NEAR(f_norm, 0.0, BASIC_TOLERANCE);
  }

TEST_F(PlasticBeamTests, RigidRotateCCWResistanceForces) {
  std::pair<real, real> deltas = rotate_ccw(in_nodes, 1.0);
  my_beam->update_state();

  vec R = my_beam->get_element_resistance_forces();
  
  EXPECT_NEAR(R(0), 0.0, BASIC_TOLERANCE); // axial
  EXPECT_NEAR(R(2), 0.0, BASIC_TOLERANCE); // vertical
  EXPECT_NEAR(R(5), 0.0, BASIC_TOLERANCE); // moment

  EXPECT_NEAR(R(6), 0.0, BASIC_TOLERANCE); // axial
  EXPECT_NEAR(R(8), 0.0, BASIC_TOLERANCE); // vertical
  EXPECT_NEAR(R(11), 0.0, BASIC_TOLERANCE); // moment
}

TEST_F(PlasticBeamTests, ConstantCompressionLength) {
  real delta = PLASTIC_BEAM_LENGTH*YIELD_STRENGTH/YOUNGS_MODULUS;
  constant_compression(in_nodes, delta);
  my_beam->update_state();

  // Calculate norms and perform assertions
  real deformed_length = my_beam->get_L();
  EXPECT_NEAR(deformed_length, PLASTIC_BEAM_LENGTH - delta, BASIC_TOLERANCE);
}

TEST_F(PlasticBeamTests, ConstantCompressionEps) {
  real delta = PLASTIC_BEAM_LENGTH*YIELD_STRENGTH/YOUNGS_MODULUS;
  constant_compression(in_nodes, delta);
  my_beam->update_state();

  // Calculate norms and perform assertions
  vec eps = my_beam->get_eps();
  EXPECT_NEAR(eps(0), -delta/PLASTIC_BEAM_LENGTH, BASIC_TOLERANCE);
  EXPECT_NEAR(eps(1), 0.0, BASIC_TOLERANCE);
}

TEST_F(PlasticBeamTests, ConstantCompressionStress) {
  real delta = PLASTIC_BEAM_LENGTH*YIELD_STRENGTH/YOUNGS_MODULUS;
  constant_compression(in_nodes, delta);
  my_beam->update_state();

  // Calculate norms and perform assertions
  vec stress = my_beam->get_local_stresses();
  EXPECT_NEAR(stress(0), -YIELD_STRENGTH*common.correct_area, BASIC_TOLERANCE);
  EXPECT_NEAR(stress(1), 0.0, BASIC_TOLERANCE);
}

TEST_F(PlasticBeamTests, ConstantCompressionLocalNodalForces) {
  real delta = PLASTIC_BEAM_LENGTH*YIELD_STRENGTH/YOUNGS_MODULUS;
  constant_compression(in_nodes, delta);
  my_beam->update_state();

  // Calculate norms and perform assertions
  vec local_f = my_beam->get_local_f();
  real local_f_norm = local_f.lpNorm<1>();
  EXPECT_NEAR(local_f(0), -YIELD_STRENGTH*common.correct_area, BASIC_TOLERANCE);
  EXPECT_NEAR(local_f(1), 0.0, BASIC_TOLERANCE);
  EXPECT_NEAR(local_f(2), 0.0, BASIC_TOLERANCE);
}

TEST_F(PlasticBeamTests, ConstantCompressionGlobalNodalForces) {
  real delta = PLASTIC_BEAM_LENGTH*YIELD_STRENGTH/YOUNGS_MODULUS;
  constant_compression(in_nodes, delta);
  my_beam->update_state();

  // Calculate norms and perform assertions
  vec R = my_beam->get_element_resistance_forces();
  real R_norm = R.lpNorm<1>();
  EXPECT_NEAR(R(0), ((delta/PLASTIC_BEAM_LENGTH)*YOUNGS_MODULUS*common.correct_area), BASIC_TOLERANCE);
  EXPECT_NEAR(R(6), -((delta/PLASTIC_BEAM_LENGTH)*YOUNGS_MODULUS*common.correct_area), BASIC_TOLERANCE);
  // the following test ensures only axial are the only forces via checking the norm.
  EXPECT_NEAR(R_norm, ((delta/PLASTIC_BEAM_LENGTH)*YOUNGS_MODULUS*common.correct_area)*2.0, BASIC_TOLERANCE); 
}

TEST_F(PlasticBeamTests, ConstantTensionLength) {
  real delta = PLASTIC_BEAM_LENGTH*YIELD_STRENGTH/YOUNGS_MODULUS;
  constant_tension(in_nodes, delta);
  my_beam->update_state();

  // Calculate norms and perform assertions
  real deformed_length = my_beam->get_L();
  EXPECT_NEAR(deformed_length, PLASTIC_BEAM_LENGTH + delta, BASIC_TOLERANCE);
}


TEST_F(PlasticBeamTests, ConstantTensionEps) {
  real delta = PLASTIC_BEAM_LENGTH*YIELD_STRENGTH/YOUNGS_MODULUS;
  constant_tension(in_nodes, delta);
  my_beam->update_state();

  // Calculate norms and perform assertions
  vec eps = my_beam->get_eps();
  EXPECT_NEAR(eps(0), (delta/PLASTIC_BEAM_LENGTH), BASIC_TOLERANCE);
  EXPECT_NEAR(eps(1), 0.0, BASIC_TOLERANCE);
}

TEST_F(PlasticBeamTests, ConstantTensionStress) {
  real delta = PLASTIC_BEAM_LENGTH*YIELD_STRENGTH/YOUNGS_MODULUS;
  constant_tension(in_nodes, delta);
  my_beam->update_state();

  // Calculate norms and perform assertions
  vec stress = my_beam->get_local_stresses();
  EXPECT_NEAR(stress(0), (delta/PLASTIC_BEAM_LENGTH)*YOUNGS_MODULUS*common.correct_area, BASIC_TOLERANCE);
  EXPECT_NEAR(stress(1), 0.0, BASIC_TOLERANCE);
}

TEST_F(PlasticBeamTests, ConstantTensionLocalNodalForces) {
  real delta = PLASTIC_BEAM_LENGTH*YIELD_STRENGTH/YOUNGS_MODULUS;
  constant_tension(in_nodes, delta);
  my_beam->update_state();

  // Calculate norms and perform assertions
  vec local_f = my_beam->get_local_f();
  real local_f_norm = local_f.lpNorm<1>();
  EXPECT_NEAR(local_f(0), YIELD_STRENGTH*common.correct_area, BASIC_TOLERANCE);
  EXPECT_NEAR(local_f(1), 0.0, BASIC_TOLERANCE);
  EXPECT_NEAR(local_f(2), 0.0, BASIC_TOLERANCE);
}

TEST_F(PlasticBeamTests, ConstantTensionGlobalNodalForces) {
  real delta = PLASTIC_BEAM_LENGTH*YIELD_STRENGTH/YOUNGS_MODULUS;
  constant_tension(in_nodes, delta);
  my_beam->update_state();

  // Calculate norms and perform assertions
  vec R = my_beam->get_element_resistance_forces();
  real R_norm = R.lpNorm<1>();
  EXPECT_NEAR(R(0), -((delta/PLASTIC_BEAM_LENGTH)*YOUNGS_MODULUS*common.correct_area), BASIC_TOLERANCE);
  EXPECT_NEAR(R(6), ((delta/PLASTIC_BEAM_LENGTH)*YOUNGS_MODULUS*common.correct_area), BASIC_TOLERANCE);
  EXPECT_NEAR(R_norm, ((delta/PLASTIC_BEAM_LENGTH)*YOUNGS_MODULUS*common.correct_area)*2.0, BASIC_TOLERANCE);
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
//   EXPECT_NEAR(eps(0), 0.0, BASIC_TOLERANCE);
//   EXPECT_NEAR(eps(1), (2/BEAM_LENGTH), BASIC_TOLERANCE);
// }

// TEST_F(PlasticBeamTests, ConstantRotationStress) {
//   constant_positive_bending(in_nodes);
//   my_beam->update_state();

//   // Calculate norms and perform assertions
//   vec stress = my_beam->get_local_stresses();
//   EXPECT_NEAR(stress(0), 0.0, BASIC_TOLERANCE);
//   EXPECT_NEAR(stress(1), (2*YOUNGS_MODULUS*SECTION_I/BEAM_LENGTH), BASIC_TOLERANCE);
// }

// TEST_F(PlasticBeamTests, ConstantRotationLocalNodalForces) {
//   constant_positive_bending(in_nodes);
//   my_beam->update_state();

//   // Calculate norms and perform assertions
//   vec local_f = my_beam->get_local_f();
//   real local_f_norm = local_f.lpNorm<1>();
//   EXPECT_NEAR(local_f(2), -(2*YOUNGS_MODULUS*SECTION_I/BEAM_LENGTH), BASIC_TOLERANCE);
//   EXPECT_NEAR(local_f(5), (2*YOUNGS_MODULUS*SECTION_I/BEAM_LENGTH), BASIC_TOLERANCE);
//   EXPECT_NEAR(local_f_norm, (2*YOUNGS_MODULUS*SECTION_I/BEAM_LENGTH)*2.0, BASIC_TOLERANCE);
// }

// TEST_F(PlasticBeamTests, ConstantRotationGlobalNodalForces) {
//   constant_positive_bending(in_nodes);
//   my_beam->update_state();

//   // Calculate norms and perform assertions
//   vec R = my_beam->get_element_resistance_forces();
//   real R_norm = R.lpNorm<1>();
//   EXPECT_NEAR(R(5), -(2*YOUNGS_MODULUS*SECTION_I/BEAM_LENGTH), BASIC_TOLERANCE);
//   EXPECT_NEAR(R(11), (2*YOUNGS_MODULUS*SECTION_I/BEAM_LENGTH), BASIC_TOLERANCE);
//   EXPECT_NEAR(R_norm, (2*YOUNGS_MODULUS*SECTION_I/BEAM_LENGTH)*2.0, BASIC_TOLERANCE);
// }

#endif