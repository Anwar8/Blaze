#ifndef TEST_HELPERS
#define TEST_HELPERS
/**
 * @name common inclusions
 * @brief each test file should only include this file. All other includes are in this file since helper functions need access to a lot of the files from Blaze.
 */
//@{      
#include <iostream>
#include <vector>
#include <numeric>
#include <cmath>
#include "gtest/gtest.h"

#include "blaze_config.hpp"
#include "ElementTypes.hpp"

#include "maths_defaults.hpp"
#include "node.hpp"

#include "BeamElementBaseClass.hpp"
#include "BeamElementCommonInterface.hpp"
#include "Linear2DBeamElement.hpp"
#include "Nonlinear2DBeamElement.hpp"


#include "ElasticPlasticMaterial.hpp"
#include "MaterialFibre.hpp"
#include "BeamColumnFiberSection.hpp"
#include "Nonlinear2DPlasticBeamElement.hpp"

#include "NodalLoad.hpp"
#include "Scribe.hpp"
#include "Model.hpp"
//@}
/**
 * @name basic definitions
 * @brief the definitions that govern all tests.
 */
//@{
// error tolerances   
#define BASIC_TOLERANCE 1e-6
#define PERCENT_TOLERANCE 0.02

// section and material properties
#define SECTION_AREA 12.437e-3
#define SECTION_I 453.266e-6
#define YIELD_STRENGTH 550e6
#define YOUNGS_MODULUS 2.06e11
#define HARDENING_RATIO_MAT 0.02
// zero hardening to allow us to calculate section properties exactly.
#define HARDENING_RATIO_FIBRE 0.0

// geometry
#define PLASTIC_BEAM_LENGTH 3.0
#define ELASTIC_BEAM_LENGTH 3.0
//@}

/**
 * @name element and section factories
 * @brief functions that create elements and sections.
 */
//@{
void initialise_I_section(BeamColumnFiberSection& section, ElasticPlasticMaterial& steel, real offset, real tf, real b, real tw, real h, int flange_divisions, int web_divisions)
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

/**
 * @brief common definitions of section parameters and many associated properties that are analytically correct.
 */
struct CommonSectionDefinitions {
    ElasticPlasticMaterial steel = ElasticPlasticMaterial(YOUNGS_MODULUS, YIELD_STRENGTH, 0.0);
    BeamColumnFiberSection I_section;
    real tf = 19.6e-3;
    real tw = 11.4e-3;
    real b = 192.8e-3;
    real h = 467.2e-3;
    real d = h - 2*tf;
    int flange_divisions = 4;
    int web_divisions = 10;
    real correct_area = tf*b*2 + (h - 2*tf)*tw; // m2
    real moment_of_inertia = tw*pow(h - 2*tf, 3)/12 + 2*b*pow(tf,3)/12 + 2*(tf*b)*pow(0.5*h - 0.5*tf, 2);
    real section_modulus = moment_of_inertia/(h/2);
    real correct_elastic_moment = section_modulus * YIELD_STRENGTH;
    real correct_plastic_moment = YIELD_STRENGTH*(tf*b)*(h - tf) + YIELD_STRENGTH*((0.5*h - tf)*tw)*(0.5*d);
    real kappa_elastic = correct_elastic_moment/(YOUNGS_MODULUS * moment_of_inertia);
    real distance_to_first_fibre = (d/40)*0.5;
    real kappa_plastic = YIELD_STRENGTH/(YOUNGS_MODULUS*distance_to_first_fibre);
    void initialise_section() 
    {
        initialise_I_section(I_section, steel, 0.0, tf, b, tw, h, 10, 40);
    }
};

void common_plastic_beam_setup(real length, std::vector<std::shared_ptr<Node>>& in_nodes, BeamColumnFiberSection& sect, std::shared_ptr<Nonlinear2DPlasticBeamElement>& my_beam, vec& U) 
{
    // Create the nodes
    in_nodes.push_back(std::make_shared<Node>(0.0, 0.0, 0.0));
    in_nodes.push_back(std::make_shared<Node>(length, 0.0, 0.0));
    my_beam = std::make_shared<Nonlinear2DPlasticBeamElement>(0, in_nodes, sect);
    // Create the d vector
    U = make_xd_vec(12);
}

void common_beam_setup(std::vector<std::shared_ptr<Node>>& in_nodes, std::shared_ptr<BeamElementBaseClass<BasicSection>>& my_beam, vec& U) 
{
    // Create the nodes
    in_nodes.push_back(std::make_shared<Node>(0.0, 0.0, 0.0));
    in_nodes.push_back(std::make_shared<Node>(ELASTIC_BEAM_LENGTH, 0.0, 0.0));
    BasicSection sect(YOUNGS_MODULUS, SECTION_AREA, SECTION_I);
    my_beam = std::make_shared<Linear2DBeamElement>(0, in_nodes, sect);
    // Create the d vector
    U = make_xd_vec(12);
}
//@}

/**
 * @name element movement
 * @brief functions that change elemental displacements.
 */
//@{
void move_laterally(std::vector<std::shared_ptr<Node>>& in_nodes, real delta = 1.0) 
{
    in_nodes[0]->set_nodal_displacement(0, delta); // node 1 U1
    in_nodes[1]->set_nodal_displacement(0, delta); // node 2 U1
}
void move_vertically(std::vector<std::shared_ptr<Node>>& in_nodes, real delta = 1.0) 
{
    in_nodes[0]->set_nodal_displacement(2, delta); // node 1 U2
    in_nodes[1]->set_nodal_displacement(2, delta); // node 2 U2
}

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

void rotate_ccw_linearly(std::vector<std::shared_ptr<Node>>& in_nodes, real theta = 2.0/ELASTIC_BEAM_LENGTH) 
{
    in_nodes[0]->set_nodal_displacement(2, -1.0); // node 1 U2
    in_nodes[0]->set_nodal_displacement(5, 2.0/ELASTIC_BEAM_LENGTH); // node 1 U33
    in_nodes[1]->set_nodal_displacement(2, 1.0); // node 2 U2
    in_nodes[1]->set_nodal_displacement(5, 2.0/ELASTIC_BEAM_LENGTH); // node 2 U33
}

void constant_compression(std::vector<std::shared_ptr<Node>>& in_nodes, real delta = 1.0) 
{
    in_nodes[0]->set_nodal_displacement(0, delta/2); // node 1 U1
    in_nodes[1]->set_nodal_displacement(0, -delta/2); // node 2 U1
}

void constant_tension(std::vector<std::shared_ptr<Node>>& in_nodes, real delta = 1.0) 
{
    in_nodes[0]->set_nodal_displacement(0, -delta/2); // node 1 U1
    in_nodes[1]->set_nodal_displacement(0, delta/2); // node 2 U1
}

void constant_positive_bending(std::vector<std::shared_ptr<Node>>& in_nodes, real theta = 1.0) 
{
    in_nodes[0]->set_nodal_displacement(5, -theta); // node 2 U33
    in_nodes[1]->set_nodal_displacement(5, theta); // node 2 U33
}
//@}


#endif