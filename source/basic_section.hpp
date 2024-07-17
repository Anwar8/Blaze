/**
 * @file basic_section.hpp
 * @brief definitions for basic sections
 * 
 */
#ifndef BASIC_SECTION_HPP
#define BASIC_SECTION_HPP

#include "maths_defaults.hpp"

/**
 * @brief basic cross-section containing geometry and material information and functions to get E, A, and I
 * 
 * @details This class currently only creates an elastic section with the equivalent properties of a 457 x 191 x 98
 * I-section made from steel. 
 * 
 */

// 457 x 191 x 98
class BasicSection {
    private:
        real E; // 2.06e11 Pa
        real A; //  0.0125 m^2
        real I; // 0.0004570000 m^4
    public:
        BasicSection() = default;
        BasicSection(real youngs_modulus, real area, real moment_of_inertia) 
        : E(youngs_modulus), A(area), I(moment_of_inertia) {}
        
        void set_E(real E) {this->E = E;}
        void set_I(real I) {this->I = I;}
        void set_A(real A) {this->A = A;}
        real get_E() {return E;}
        real get_A() {return A;}
        real get_I() {return I;}
};
#endif