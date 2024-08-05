#ifndef ELASTIC_PLASTIC_MATERIAL_TESTS_HPP
#define ELASTIC_PLASTIC_MATERIAL_TESTS_HPP

#include "gtest/gtest.h"
#include "../ElasticPlasticMaterial.hpp"
#define YOUNGS_MODULUS_MAT 2e11
#define YIELD_STRENGTH_MAT 355e6
#define HARDENING_RATIO 0.02
#define TOLERANCE 1e-6


class ElasticPlasticMaterialTest : public ::testing::Test {
  public:
    ElasticPlasticMaterial Steel = ElasticPlasticMaterial(YOUNGS_MODULUS_MAT, YIELD_STRENGTH_MAT, HARDENING_RATIO*YOUNGS_MODULUS_MAT);
    real yield_strain = YIELD_STRENGTH_MAT/YOUNGS_MODULUS_MAT;
  
    void SetUp() override {
        
    }
    void TearDown() override {
}
};

/**
 * @details checks that \ref ElasticPlasticMaterial correctly identifies still elastic state and stress associated with it in tension.
 * 
 */
TEST_F(ElasticPlasticMaterialTest, IncrementStrainElastic)
{
    Steel.increment_strain(0.9*yield_strain);
    real correct_stress = YOUNGS_MODULUS_MAT*0.9*yield_strain;
    real calculated_stress = Steel.get_stress();
    EXPECT_TRUE(Steel.is_elastic());
    EXPECT_NEAR(calculated_stress, correct_stress, TOLERANCE);
}

/**
 * @details checks that \ref ElasticPlasticMaterial correctly identifies still elastic state and stress associated with it in compression.
 * 
 */
TEST_F(ElasticPlasticMaterialTest, IncrementStrainElasticCompression)
{
    Steel.increment_strain(-0.9*yield_strain);
    real correct_stress = -YOUNGS_MODULUS_MAT*0.9*yield_strain;
    real calculated_stress = Steel.get_stress();
    EXPECT_TRUE(Steel.is_elastic());
    EXPECT_NEAR(calculated_stress, correct_stress, TOLERANCE);
}

/**
 * @details checks that \ref ElasticPlasticMaterial correctly identifies yielded state and stress associated with it in tension.
 * 
 */
TEST_F(ElasticPlasticMaterialTest, IncrementStrainPlastic)
{
    Steel.increment_strain(1.1*yield_strain);
    real E_tangent =  YOUNGS_MODULUS_MAT * HARDENING_RATIO*YOUNGS_MODULUS_MAT / (YOUNGS_MODULUS_MAT + HARDENING_RATIO*YOUNGS_MODULUS_MAT);
    real correct_stress = YIELD_STRENGTH_MAT + E_tangent*0.1*yield_strain;

    real calculated_stress = Steel.get_stress();
    EXPECT_FALSE(Steel.is_elastic());
    EXPECT_NEAR(calculated_stress, correct_stress, TOLERANCE);
}

/**
 * @details verifies that the \ref ElasticPlasticMaterial is able to correctly unload after loading post-yielding in tension, and then
 * unloading in compression.
 */
TEST_F(ElasticPlasticMaterialTest, IncrementStrainElasticUnloading)
{
    Steel.increment_strain(1.1*yield_strain);
    real E_tangent =  YOUNGS_MODULUS_MAT * HARDENING_RATIO*YOUNGS_MODULUS_MAT / (YOUNGS_MODULUS_MAT + HARDENING_RATIO*YOUNGS_MODULUS_MAT);
    real max_stress = YIELD_STRENGTH_MAT + E_tangent*0.1*yield_strain;
    Steel.update_starting_state();

    Steel.increment_strain(-0.2*yield_strain);
    real correct_stress = max_stress -0.2*yield_strain*YOUNGS_MODULUS_MAT;
    real calculated_stress = Steel.get_stress();
    EXPECT_TRUE(Steel.is_elastic());
    EXPECT_NEAR(calculated_stress, correct_stress, TOLERANCE);
}

/**
 * @details verifies that the \ref ElasticPlasticMaterial object correctly updates its internal state when subjected to a tensile strain in the plastic range, 
 * and that the calculated plastic strain and yield stress match the expected values.
 * 
 */
TEST_F(ElasticPlasticMaterialTest, IncrementStrainPlasticStrainYieldStrength)
{
    Steel.increment_strain(1.1*yield_strain);
    real E_tangent =  YOUNGS_MODULUS_MAT * HARDENING_RATIO*YOUNGS_MODULUS_MAT / (YOUNGS_MODULUS_MAT + HARDENING_RATIO*YOUNGS_MODULUS_MAT);
    real max_stress = YIELD_STRENGTH_MAT + E_tangent*0.1*yield_strain;
    
    Steel.evolve_yield_surface();
    

    real beta = (yield_strain + (max_stress - YIELD_STRENGTH_MAT)/YOUNGS_MODULUS_MAT)/(1.1*yield_strain);
    real correct_plastic_strain = (1.1*yield_strain)*(1 - beta);
    real correct_fy_bar =  YIELD_STRENGTH_MAT + correct_plastic_strain*HARDENING_RATIO*YOUNGS_MODULUS_MAT;

    real calculated_plastic_strain = Steel.get_plastic_strain();
    real calculated_fy_bar = Steel.get_fy_bar();

    EXPECT_NEAR(calculated_plastic_strain, correct_plastic_strain, TOLERANCE);
    EXPECT_NEAR(calculated_fy_bar, correct_fy_bar, TOLERANCE);
}

/**
 * @details checks that \ref ElasticPlasticMaterial correctly identifies yielded state and stress associated with it in compression.
 * 
 */
TEST_F(ElasticPlasticMaterialTest, IncrementStrainPlasticCompression)
{
    Steel.increment_strain(-1.1*yield_strain);
    real E_tangent =  YOUNGS_MODULUS_MAT * HARDENING_RATIO*YOUNGS_MODULUS_MAT / (YOUNGS_MODULUS_MAT + HARDENING_RATIO*YOUNGS_MODULUS_MAT);
    real correct_stress = -YIELD_STRENGTH_MAT - E_tangent*0.1*yield_strain;

    real calculated_stress = Steel.get_stress();
    EXPECT_FALSE(Steel.is_elastic());
    EXPECT_NEAR(calculated_stress, correct_stress, TOLERANCE);
}

/**
 * @details verifies that the \ref ElasticPlasticMaterial is able to correctly unload after loading post-yielding in compression, and then
 * unloading in tension.
 */
TEST_F(ElasticPlasticMaterialTest, IncrementStrainElasticUnloadingCompression)
{
    Steel.increment_strain(-1.1*yield_strain);
    real E_tangent =  YOUNGS_MODULUS_MAT * HARDENING_RATIO*YOUNGS_MODULUS_MAT / (YOUNGS_MODULUS_MAT + HARDENING_RATIO*YOUNGS_MODULUS_MAT);
    real max_stress = -YIELD_STRENGTH_MAT - E_tangent*0.1*yield_strain;
    Steel.update_starting_state();

    Steel.increment_strain(0.2*yield_strain);
    real correct_stress = max_stress +0.2*yield_strain*YOUNGS_MODULUS_MAT;
    real calculated_stress = Steel.get_stress();
    EXPECT_TRUE(Steel.is_elastic());
    EXPECT_NEAR(calculated_stress, correct_stress, TOLERANCE);
}

/**
 * @details verifies that the \ref ElasticPlasticMaterial object correctly updates its internal state when subjected to a compressive strain in the plastic range, 
 * and that the calculated plastic strain and yield stress match the expected values.
 * 
 */
TEST_F(ElasticPlasticMaterialTest, IncrementStrainPlasticStrainYieldStrengthCompression)
{
    Steel.increment_strain(-1.1*yield_strain);
    real E_tangent =  YOUNGS_MODULUS_MAT * HARDENING_RATIO*YOUNGS_MODULUS_MAT / (YOUNGS_MODULUS_MAT + HARDENING_RATIO*YOUNGS_MODULUS_MAT);
    real max_stress = -YIELD_STRENGTH_MAT - E_tangent*0.1*yield_strain;
    
    Steel.evolve_yield_surface();
    

    real beta = (yield_strain + (abs(max_stress) - YIELD_STRENGTH_MAT)/YOUNGS_MODULUS_MAT)/(1.1*yield_strain);
    real correct_plastic_strain = (1.1*yield_strain)*(1 - beta);
    real correct_fy_bar =  YIELD_STRENGTH_MAT + correct_plastic_strain*HARDENING_RATIO*YOUNGS_MODULUS_MAT;

    real calculated_plastic_strain = Steel.get_plastic_strain();
    real calculated_fy_bar = Steel.get_fy_bar();

    EXPECT_NEAR(calculated_plastic_strain, correct_plastic_strain, TOLERANCE);
    EXPECT_NEAR(calculated_fy_bar, correct_fy_bar, TOLERANCE);
}

/**
 * @brief Checks that loading and unloading all strain results in 0.0 strain for the material, but nonzero stress for the material.
 * 
 */
TEST_F(ElasticPlasticMaterialTest, IncrementStrainCyclicZeroStrain)
{
    Steel.increment_strain(1.1*yield_strain);
    Steel.update_starting_state();
    Steel.increment_strain(-1.1*yield_strain);
    
    real calculated_stress = Steel.get_stress();
    real calculated_strain = Steel.get_strain();

    EXPECT_TRUE(Steel.is_elastic());
    EXPECT_NEAR(calculated_strain, 0.0, TOLERANCE);
    EXPECT_TRUE(calculated_stress != 0.0);
}

/**
 * @brief Checks that loading and unloading ALL elastic strain results in 0.0 stress for the material, and positive strains.
 * 
 */
TEST_F(ElasticPlasticMaterialTest, IncrementStrainCyclicZeroStress)
{
    Steel.increment_strain(1.1*yield_strain);
    Steel.update_starting_state();

    Steel.evolve_yield_surface();
    real new_yield_strain = Steel.get_fy_bar()/YOUNGS_MODULUS_MAT;
    Steel.increment_strain(-new_yield_strain);
    
    real calculated_stress = Steel.get_stress();
    real calculated_strain = Steel.get_strain();

    EXPECT_TRUE(Steel.is_elastic());
    EXPECT_NEAR(calculated_stress, 0.0, TOLERANCE);
    EXPECT_TRUE(calculated_strain > 0);
    
}


/**
 * @brief Checks that loading and unloading all strain results in non-zero accumulated plastic strain
 * as isotropic hardening links hardening to accumulated plastic strain.
 * 
 */
TEST_F(ElasticPlasticMaterialTest, IncrementStrainCyclicNonZeroPlastic)
{
    Steel.increment_strain(1.1*yield_strain);
    real E_tangent =  YOUNGS_MODULUS_MAT * HARDENING_RATIO*YOUNGS_MODULUS_MAT / (YOUNGS_MODULUS_MAT + HARDENING_RATIO*YOUNGS_MODULUS_MAT);
    real max_stress = YIELD_STRENGTH_MAT + E_tangent*0.1*yield_strain;
    real beta = (yield_strain + (abs(max_stress) - YIELD_STRENGTH_MAT)/YOUNGS_MODULUS_MAT)/(1.1*yield_strain);
    real correct_plastic_strain = (1.1*yield_strain)*(1 - beta);
    Steel.update_starting_state();
    Steel.increment_strain(-1.1*yield_strain);

    real calculated_plastic_strain = Steel.get_plastic_strain();
    EXPECT_NEAR(calculated_plastic_strain, correct_plastic_strain, TOLERANCE);
}

#endif