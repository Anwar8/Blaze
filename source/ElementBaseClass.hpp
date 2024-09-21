/**
 * @file ElementBaseClass.hpp
 * @brief empty class the serves as a placeholder in containers.
 * 
 */
#ifndef ELEMENT_BASE_CLASS_HPP
#define ELEMENT_BASE_CLASS_HPP
#include "maths_defaults.hpp"
class ElementBaseClass 
{
    public:
        virtual void map_stiffness() = 0;
        virtual void calc_global_stiffness_triplets() = 0;
        virtual void update_state() = 0;
        virtual void print_info() = 0;
        virtual void print_element_state(bool print_stresses = true, bool print_strains = false,
                                 bool print_nodal_disp = false, bool print_nodal_forces = false) = 0;
        virtual std::vector<spnz> get_global_resistance_force_triplets() = 0;
        virtual std::vector<spnz> get_global_stiffness_triplets() = 0;

};
#endif