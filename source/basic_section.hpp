/**
 * @file basic_section.hpp
 * @brief definitions for basic sections
 * 
 */
#ifndef BASIC_SECTION_HPP
#define BASIC_SECTION_HPP

#include "maths_defaults.hpp"
#include "SectionBaseClass.hpp"

/**
 * @brief basic cross-section containing geometry and material information and functions to get E, A, and I
 * 
 */

class BasicSection : public SectionBaseClass {
    private:
        real E;
        real A;
        real I;
    public:
        BasicSection() = default;
        BasicSection(real youngs_modulus, real area, real moment_of_inertia) 
        : E(youngs_modulus), A(area), I(moment_of_inertia) {section_type = Basic;}
        
        void set_E(real E) {this->E = E;}
        void set_I(real I) {this->I = I;}
        void set_A(real A) {this->A = A;}
        real get_E() {return E;}
        real get_A() {return A;}
        real get_I() {return I;}
};
#endif