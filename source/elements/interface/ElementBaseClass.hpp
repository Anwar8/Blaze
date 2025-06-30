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
        virtual void update_section_starting_state() = 0;
        virtual void print_info() = 0;
        virtual void print_element_state(bool print_stresses = true, bool print_strains = false,
                                 bool print_nodal_disp = false, bool print_nodal_forces = false) = 0;
        virtual std::vector<spnz> get_global_resistance_force_triplets() = 0;
        /**
         * @brief inserts the contents of \ref element_global_resistance_forces into the end of global resistance triplets vector \ref Assembler::R_global_triplets. 
         * Used to reduce copying during assembly, still basically a getter function.
         * 
         * @param R_global_triplets The container for the triplets that are used for assembling the global resistance vector \f$ \boldsymbol{R}\f$ from element contributions.
         */ 
        virtual void insert_global_resistance_force_triplets(std::vector<spnz>& global_resistance_triplets_vector) = 0;
        virtual std::vector<spnz> get_global_stiffness_triplets() = 0;
        /**
         * @brief inserts the contents of \ref global_stiffness_triplets into the end of global_triplets_vector. Used to reduce copying during assembly, still basically a getter function.
         * 
         * @param global_triplets_vector The container for the triplets that are used for assembling the global stiffness matrix \f$ \boldsymbol{K}\f$ from element contributions.
         */
        virtual void insert_global_stiffness_triplets(std::vector<spnz>& global_triplets_vector) = 0;
        virtual size_t get_id() const = 0;
        

};
#endif