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
        BasicSection() {
            section_type = Basic;
        };
        BasicSection(real youngs_modulus, real area, real moment_of_inertia) 
        : E(youngs_modulus), A(area), I(moment_of_inertia) { section_type = Basic;}

        BasicSection(const BasicSection& other)
        {
            E = other.E;
            A = other.A;
            I = other.I;
            section_type = other.section_type;
        }
        
        void set_E(real E) {this->E = E;}
        void set_I(real I) {this->I = I;}
        void set_A(real A) {this->A = A;}
        
        real get_E() {return E;}
        real get_A() {return A;}
        real get_I() {return I;}

        virtual void update_section_state(vec& eps) {};
        mat get_D_t() {return make_xd_mat(2,2);}
        
};
#endif