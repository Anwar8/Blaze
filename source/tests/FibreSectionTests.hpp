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

// class FibreSectionPureBendingTests : public ::testing::Test {
//   public:
//     ElasticPlasticMaterial steel = ElasticPlasticMaterial(YOUNGS_MODULUS_FIBRE, YIELD_STRENGTH_FIBRE, HARDENING_RATIO_FIBRE*YOUNGS_MODULUS_FIBRE);
//     real yield_strain = YIELD_STRENGTH_FIBRE/YOUNGS_MODULUS_FIBRE;
  
//     void SetUp() override {
        
//     }
//     void TearDown() override {
// }
// };

#endif 