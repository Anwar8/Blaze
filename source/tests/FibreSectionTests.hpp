#ifndef FIBRE_SECTIONTESTS_HPP
#define FIBRE_SECTIONTESTS_HPP

#include "gtest/gtest.h"
#include "../ElasticPlasticMaterial.hpp"
#include "../BeamColumnFiberSection.hpp"
#include "../MaterialFibre.hpp"

#define YOUNGS_MODULUS_FIBRE 2e11
#define YIELD_STRENGTH_FIBRE 355e6
// zero hardening to allow us to calculate section properties exactly.
#define HARDENING_RATIO_FIBRE 0.0
#define TOLERANCE_FIBRE 1e-6
#define TOLERANCE_SECTION_FORCES 0.02


void build_I_section(BeamColumnFiberSection& section, ElasticPlasticMaterial& steel, real offset, real tf, real b, real tw, real h, int flange_divisions, int web_divisions)
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

class MaterialFibreTests : public ::testing::Test {
  public:
    ElasticPlasticMaterial steel = ElasticPlasticMaterial(YOUNGS_MODULUS_FIBRE, YIELD_STRENGTH_FIBRE, HARDENING_RATIO_FIBRE*YOUNGS_MODULUS_FIBRE);
    MaterialFibre fibre;
    void SetUp() override {
        fibre.initialise_fibre(&steel, 1.0, 0.5, 0.0);
    }
    void TearDown() override {

}
};

TEST_F(MaterialFibreTests, AreaCorrect)
{
    real saved_area = fibre.get_area();
    EXPECT_NEAR(saved_area, 1.0, TOLERANCE_FIBRE);
}

class FibreSectionCentroidTests : public ::testing::Test {
  public:
    ElasticPlasticMaterial steel = ElasticPlasticMaterial(YOUNGS_MODULUS_FIBRE, YIELD_STRENGTH_FIBRE, HARDENING_RATIO_FIBRE*YOUNGS_MODULUS_FIBRE);
    BeamColumnFiberSection I_section;
    // UB457 x 191 x 98: h: 467.2	b:192.8	tw:11.4	tf:19.6
    
    real tf = 19.6e-3;
    real tw = 11.4e-3;
    real b = 192.8e-3;
    real h = 467.2e-3;
    real correct_area = tf*b*2 + (h - 2*tf)*tw; // m2
    void SetUp() override {
        
    }
    void TearDown() override {
}
};

/**
 * @brief checks that \ref BeamColumnFiberSection correctly calculates the area of the section.
 * 
 */
TEST_F(FibreSectionCentroidTests, CalculatedArea)
{
    build_I_section(I_section, steel, 0.0, tf, b, tw, h, 10, 40);
    I_section.calc_area_weighted_E();
    I_section.calc_section_centroid();
    real calculated_area = I_section.get_section_area();
    
    EXPECT_NEAR(calculated_area, correct_area, TOLERANCE_FIBRE);
}

/**
 * @brief checks that \ref BeamColumnFiberSection correctly calculates the centroid of the section.
 * 
 */
TEST_F(FibreSectionCentroidTests, CalculatedCentroidNoOffset)
{
    build_I_section(I_section, steel, 0.0, tf, b, tw, h, 10, 40);
    I_section.calc_area_weighted_E();
    I_section.calc_section_centroid();
    real calculated_section_centroid = I_section.get_y_bar();
    real correct_section_centroid = h/2;
    
    EXPECT_NEAR(calculated_section_centroid, correct_section_centroid, TOLERANCE_FIBRE);
}


/**
 * @brief checks that \ref BeamColumnFiberSection correctly calculates the centroid of the section if there is a positive offset.
 * 
 */
TEST_F(FibreSectionCentroidTests, CalculatedCentroidPositiveOffset)
{
    real offset = 2.3;
    build_I_section(I_section, steel, offset, tf, b, tw, h, 10, 40);
    I_section.calc_area_weighted_E();
    I_section.calc_section_centroid();
    real calculated_section_centroid = I_section.get_y_bar();
    real correct_section_centroid = offset + h/2;
    
    EXPECT_NEAR(calculated_section_centroid, correct_section_centroid, TOLERANCE_FIBRE);
}

/**
 * @brief checks that \ref BeamColumnFiberSection correctly calculates the centroid of the section if there is a negative offset.
 * 
 */
TEST_F(FibreSectionCentroidTests, CalculatedCentroidNegativeOffset)
{
    real offset = -1.00;
    build_I_section(I_section, steel, offset, tf, b, tw, h, 10, 40);
    I_section.calc_area_weighted_E();
    I_section.calc_section_centroid();
    real calculated_section_centroid = I_section.get_y_bar();
    real correct_section_centroid = offset + h/2;
    
    EXPECT_NEAR(calculated_section_centroid, correct_section_centroid, TOLERANCE_FIBRE);
}

/**
 * @brief checks that \ref BeamColumnFiberSection correctly calculates the centroid of the section if the real centroid is at zero.
 * 
 */
TEST_F(FibreSectionCentroidTests, CalculatedZeroCentroid)
{
    real offset = -h/2;
    build_I_section(I_section, steel, offset, tf, b, tw, h, 10, 40);
    I_section.calc_area_weighted_E();
    I_section.calc_section_centroid();
    real calculated_section_centroid = I_section.get_y_bar();
    real correct_section_centroid = offset + h/2;
    
    EXPECT_NEAR(calculated_section_centroid, correct_section_centroid, TOLERANCE_FIBRE);
}

class FibreSectionPureBendingTests : public ::testing::Test {
  public:
    ElasticPlasticMaterial steel = ElasticPlasticMaterial(YOUNGS_MODULUS_FIBRE, YIELD_STRENGTH_FIBRE, HARDENING_RATIO_FIBRE*YOUNGS_MODULUS_FIBRE);
    BeamColumnFiberSection I_section;
    
    real tf = 19.6e-3;
    real tw = 11.4e-3;
    real b = 192.8e-3;
    real h = 467.2e-3;
    real d = h - 2*tf;

    real moment_of_inertia = tw*pow(h - 2*tf, 3)/12 + 2*b*pow(tf,3)/12 + 2*(tf*b)*pow(0.5*h - 0.5*tf, 2);
    real section_modulus = moment_of_inertia/(h/2);
    real correct_elastic_moment = section_modulus * YIELD_STRENGTH_FIBRE;
    real correct_plastic_moment = YIELD_STRENGTH_FIBRE*(tf*b)*(h - tf) + YIELD_STRENGTH_FIBRE*((0.5*h - tf)*tw)*(0.5*d);
    real kappa_elastic = correct_elastic_moment/(YOUNGS_MODULUS_FIBRE * moment_of_inertia);

    real distance_to_first_fibre = (d/40)*0.5;
    real kappa_plastic = YIELD_STRENGTH_FIBRE/(YOUNGS_MODULUS_FIBRE*distance_to_first_fibre);
    void SetUp() override {
        build_I_section(I_section, steel, 0.0, tf, b, tw, h, 10, 40);
        I_section.calc_area_weighted_E();
        I_section.calc_section_centroid();
    }
    void TearDown() override {
}
};

/**
 * @brief checks that \ref BeamColumnFiberSection correctly calculates the elastic moment of a section.
 * 
 */
TEST_F(FibreSectionPureBendingTests, ElasticMoment)
{
    I_section.increment_section_strains(0.0, kappa_elastic);
    I_section.increment_fibre_strains();
    I_section.calc_section_forces();
    real calculated_section_moment = I_section.get_moment_yy();
    real calculated_section_force = I_section.get_axial_force();

    
    real percent_error = abs(correct_elastic_moment - calculated_section_moment)/abs(correct_elastic_moment);

    EXPECT_NEAR(percent_error, 0.0, TOLERANCE_SECTION_FORCES);
    EXPECT_NEAR(calculated_section_force, 0.0, TOLERANCE_FIBRE);
}

/**
 * @brief checks that \ref BeamColumnFiberSection correctly calculates the plastic moment of a section.
 * 
 */
TEST_F(FibreSectionPureBendingTests, PlasticMoment)
{
    I_section.increment_section_strains(0.0, kappa_plastic);
    I_section.increment_fibre_strains();
    I_section.calc_section_forces();
    real calculated_section_moment = I_section.get_moment_yy();
    real calculated_section_force = I_section.get_axial_force();
    
    real percent_error = abs(correct_plastic_moment - calculated_section_moment)/abs(correct_plastic_moment);

    EXPECT_NEAR(percent_error, 0.0, TOLERANCE_SECTION_FORCES);
    EXPECT_NEAR(calculated_section_force, 0.0, TOLERANCE_FIBRE);
}



/**
 * @brief checks that \ref BeamColumnFiberSection correctly calculates the negative elastic moment of a section.
 * 
 */
TEST_F(FibreSectionPureBendingTests, NegativeElasticMoment)
{
    I_section.increment_section_strains(0.0, -kappa_elastic);
    I_section.increment_fibre_strains();
    I_section.calc_section_forces();
    real calculated_section_moment = I_section.get_moment_yy();
    real calculated_section_force = I_section.get_axial_force();

    
    real percent_error = abs(-correct_elastic_moment - calculated_section_moment)/abs(-correct_elastic_moment);

    EXPECT_NEAR(percent_error, 0.0, TOLERANCE_SECTION_FORCES);
    EXPECT_NEAR(calculated_section_force, 0.0, TOLERANCE_FIBRE);
}

/**
 * @brief checks that \ref BeamColumnFiberSection correctly calculates the negative plastic moment of a section.
 * 
 */
TEST_F(FibreSectionPureBendingTests, NegativePlasticMoment)
{
    I_section.increment_section_strains(0.0, -kappa_plastic);
    I_section.increment_fibre_strains();
    I_section.calc_section_forces();
    real calculated_section_moment = I_section.get_moment_yy();
    real calculated_section_force = I_section.get_axial_force();
    
    real percent_error = abs(-correct_plastic_moment - calculated_section_moment)/abs(-correct_plastic_moment);

    EXPECT_NEAR(percent_error, 0.0, TOLERANCE_SECTION_FORCES);
    EXPECT_NEAR(calculated_section_force, 0.0, TOLERANCE_FIBRE);
}




class FibreSectionPureAxialTests : public ::testing::Test {
  public:
    ElasticPlasticMaterial steel = ElasticPlasticMaterial(YOUNGS_MODULUS_FIBRE, YIELD_STRENGTH_FIBRE, HARDENING_RATIO_FIBRE*YOUNGS_MODULUS_FIBRE);
    BeamColumnFiberSection I_section;
    
    real tf = 19.6e-3;
    real tw = 11.4e-3;
    real b = 192.8e-3;
    real h = 467.2e-3;
    real d = h - 2*tf;

    real section_area = d*tw + 2*tf*b;
    real axial_yield_strain = YIELD_STRENGTH_FIBRE/YOUNGS_MODULUS_FIBRE;
    real correct_force = YIELD_STRENGTH_FIBRE*section_area;

    void SetUp() override {
        build_I_section(I_section, steel, 0.0, tf, b, tw, h, 10, 40);
        I_section.calc_area_weighted_E();
        I_section.calc_section_centroid();
    }
    void TearDown() override {
}
};

/**
 * @brief checks that \ref BeamColumnFiberSection correctly calculates the tensile axial force of a section.
 * 
 */
TEST_F(FibreSectionPureAxialTests, YieldTensileForce)
{
    I_section.increment_section_strains(axial_yield_strain, 0.0);
    I_section.increment_fibre_strains();
    I_section.calc_section_forces();
    real calculated_section_force = I_section.get_axial_force();
    real calculated_section_moment = I_section.get_moment_yy();

    real percent_error = abs(correct_force - calculated_section_force)/abs(correct_force);

    EXPECT_NEAR(percent_error, 0.0, TOLERANCE_SECTION_FORCES);
    EXPECT_NEAR(calculated_section_moment, 0.0, TOLERANCE_FIBRE);
}

/**
 * @brief checks that \ref BeamColumnFiberSection correctly calculates the compressive axial force of a section.
 * 
 */
TEST_F(FibreSectionPureAxialTests, YieldCompressiveForce)
{
    I_section.increment_section_strains(-axial_yield_strain, 0.0);
    I_section.increment_fibre_strains();
    I_section.calc_section_forces();
    real calculated_section_force = I_section.get_axial_force();
    real calculated_section_moment = I_section.get_moment_yy();

    real percent_error = abs(-correct_force - calculated_section_force)/abs(-correct_force);

    EXPECT_NEAR(percent_error, 0.0, TOLERANCE_SECTION_FORCES);
    EXPECT_NEAR(calculated_section_moment, 0.0, TOLERANCE_FIBRE);
}


/**
 * @brief checks that \ref BeamColumnFiberSection correctly calculates the tensile axial force of a section beyond yield.
 * 
 */
TEST_F(FibreSectionPureAxialTests, PostYieldTensileForce)
{
    I_section.increment_section_strains(1.2*axial_yield_strain, 0.0);
    I_section.increment_fibre_strains();
    I_section.calc_section_forces();
    real calculated_section_force = I_section.get_axial_force();
    real calculated_section_moment = I_section.get_moment_yy();

    real percent_error = abs(correct_force - calculated_section_force)/abs(correct_force);

    EXPECT_NEAR(percent_error, 0.0, TOLERANCE_SECTION_FORCES);
    EXPECT_NEAR(calculated_section_moment, 0.0, TOLERANCE_FIBRE);
}

/**
 * @brief checks that \ref BeamColumnFiberSection correctly calculates the compressive axial force of a section beyond yield.
 * 
 */
TEST_F(FibreSectionPureAxialTests, PostYieldCompressiveForce)
{
    I_section.increment_section_strains(-1.2*axial_yield_strain, 0.0);
    I_section.increment_fibre_strains();
    I_section.calc_section_forces();
    real calculated_section_force = I_section.get_axial_force();
    real calculated_section_moment = I_section.get_moment_yy();

    real percent_error = abs(-correct_force - calculated_section_force)/abs(-correct_force);

    EXPECT_NEAR(percent_error, 0.0, TOLERANCE_SECTION_FORCES);
    EXPECT_NEAR(calculated_section_moment, 0.0, TOLERANCE_FIBRE);
}

/**
 * @brief checks that \ref BeamColumnFiberSection correctly calculates the compressive axial force of a section as being elastic if the starting state is not updated.
 * 
 */
TEST_F(FibreSectionPureAxialTests, IncrementalNoYieldCompressiveForce)
{
    vec d_eps = make_xd_vec(2);
    d_eps(0) = -0.6*axial_yield_strain;
    I_section.update_section_state(d_eps);
    I_section.update_section_state(d_eps);
    real calculated_section_force = I_section.get_axial_force();
    real calculated_section_moment = I_section.get_moment_yy();
    real non_yield_correct_force =0.6*correct_force;

    real percent_error = abs(-non_yield_correct_force - calculated_section_force)/abs(-non_yield_correct_force);

    EXPECT_NEAR(percent_error, 0.0, TOLERANCE_SECTION_FORCES);
    EXPECT_NEAR(calculated_section_moment, 0.0, TOLERANCE_FIBRE);
}

/**
 * @brief checks that \ref BeamColumnFiberSection correctly calculates the compressive axial force of a section beyond yield if applied in increments.
 * 
 */
TEST_F(FibreSectionPureAxialTests, IncrementalPostYieldCompressiveForce)
{
    vec d_eps = make_xd_vec(2);
    d_eps(0) = -0.6*axial_yield_strain;
    I_section.update_section_state(d_eps);
    I_section.update_section_starting_state();
    d_eps(0) = -1.2*axial_yield_strain;
    I_section.update_section_state(d_eps);
    // I_section.increment_section_strains(-1.2*axial_yield_strain, 0.0);
    // I_section.increment_fibre_strains();
    // I_section.calc_section_forces();
    real calculated_section_force = I_section.get_axial_force();
    real calculated_section_moment = I_section.get_moment_yy();

    real percent_error = abs(-correct_force - calculated_section_force)/abs(-correct_force);

    EXPECT_NEAR(percent_error, 0.0, TOLERANCE_SECTION_FORCES);
    EXPECT_NEAR(calculated_section_moment, 0.0, TOLERANCE_FIBRE);
}


class FibreSectionTangentConstitutiveMatrix : public ::testing::Test {
  public:
    
    ElasticPlasticMaterial steel = ElasticPlasticMaterial(YOUNGS_MODULUS_FIBRE, YIELD_STRENGTH_FIBRE, HARDENING_RATIO_FIBRE*YOUNGS_MODULUS_FIBRE);
    BeamColumnFiberSection I_section;
    
    real tf = 19.6e-3;
    real tw = 11.4e-3;
    real b = 192.8e-3;
    real h = 467.2e-3;
    real d = h - 2*tf;

    real section_area = d*tw + 2*tf*b;
    real moment_of_inertia = tw*pow(h - 2*tf, 3)/12 + 2*b*pow(tf,3)/12 + 2*(tf*b)*pow(0.5*h - 0.5*tf, 2);

    real axial_yield_strain = YIELD_STRENGTH_FIBRE/YOUNGS_MODULUS_FIBRE;
    real correct_EI = YOUNGS_MODULUS_FIBRE*moment_of_inertia;
    real correct_EA = YOUNGS_MODULUS_FIBRE*section_area;
    

    void SetUp() override {
        build_I_section(I_section, steel, 0.0, tf, b, tw, h, 10, 40);
    }
    void TearDown() override {
}
};

/**
 * @brief checks that \ref BeamColumnFiberSection correctly calculates the tangent constitutive matrix if the section is elastic.
 * 
 */
TEST_F(FibreSectionTangentConstitutiveMatrix, CheckDtNoYield)
{
    mat D_t_calculated = make_xd_mat(2,2);
    vec d_eps = make_xd_vec(2);
    
    d_eps(0) = 0.9*axial_yield_strain;
    I_section.update_section_state(d_eps);
    D_t_calculated = I_section.get_D_t();
    
    real calculated_EA = D_t_calculated(0,0);
    real calculated_EI = D_t_calculated(1,1);

    
    real percent_error_00 = abs(correct_EA - calculated_EA)/abs(correct_EA);
    real percent_error_11 = abs(correct_EI - calculated_EI)/abs(correct_EI);

    EXPECT_NEAR(percent_error_00, 0.0, TOLERANCE_SECTION_FORCES);
    EXPECT_NEAR(percent_error_11, 0.0, TOLERANCE_SECTION_FORCES);
}

class FibreSectionTangentConstitutiveMatrixHardening : public ::testing::Test {
  public:
    real hardening_param_H = 0.01*YOUNGS_MODULUS_FIBRE;
    ElasticPlasticMaterial hardening_steel = ElasticPlasticMaterial(YOUNGS_MODULUS_FIBRE, YIELD_STRENGTH_FIBRE, hardening_param_H);
    BeamColumnFiberSection I_section;
    real correct_E_t = YOUNGS_MODULUS_FIBRE * hardening_param_H/ (YOUNGS_MODULUS_FIBRE + hardening_param_H);
    real tf = 19.6e-3;
    real tw = 11.4e-3;
    real b = 192.8e-3;
    real h = 467.2e-3;
    real d = h - 2*tf;

    real section_area = d*tw + 2*tf*b;
    real moment_of_inertia = tw*pow(h - 2*tf, 3)/12 + 2*b*pow(tf,3)/12 + 2*(tf*b)*pow(0.5*h - 0.5*tf, 2);

    real axial_yield_strain = YIELD_STRENGTH_FIBRE/YOUNGS_MODULUS_FIBRE;
    real correct_EI = correct_E_t*moment_of_inertia;
    real correct_EA = correct_E_t*section_area;
    

    void SetUp() override {
        build_I_section(I_section, hardening_steel, 0.0, tf, b, tw, h, 10, 40);
    }
    void TearDown() override {
}
};

/**
 * @brief checks that \ref BeamColumnFiberSection correctly calculates the tangent constitutive matrix if the section is plastic and with hardening.
 * 
 */
TEST_F(FibreSectionTangentConstitutiveMatrixHardening, CheckDtPostYield)
{
    mat D_t_calculated = make_xd_mat(2,2);
    vec d_eps = make_xd_vec(2);
    
    d_eps(0) = 1.1*axial_yield_strain;
    I_section.update_section_state(d_eps);
    D_t_calculated = I_section.get_D_t();
    
    real calculated_EA = D_t_calculated(0,0);
    real calculated_EI = D_t_calculated(1,1);

    
    real percent_error_00 = abs(correct_EA - calculated_EA)/abs(correct_EA);
    real percent_error_11 = abs(correct_EI - calculated_EI)/abs(correct_EI);

    EXPECT_NEAR(percent_error_00, 0.0, TOLERANCE_SECTION_FORCES);
    EXPECT_NEAR(percent_error_11, 0.0, TOLERANCE_SECTION_FORCES);
}
#endif 
