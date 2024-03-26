/**
 * @file basic_shape_function.hpp
 * @brief definitions for basic shape functions and related local matrices
 * 
 */
#ifndef BASIC_SHAPE_FUNCTION
#define BASIC_SHAPE_FUNCTION

#include "maths_defaults.hpp"
#include "basic_section.hpp"

/**
 * @brief shape function class that contains the local stiffness, freedoms, derivative, and order of DoFs
 * 
 * @details the shape function class is used to calculate the values of the shape function related matrices
 * for each row columns of the N and B. This requires knowing where along the element these values will be
 * calculated and how long the element is. Given section properties, this class will also calculate local
 * stiffness. Currently, this is done assuming an elastic non-changing local stiffness matrix.
 * 
 */
class BasicShapeFunction {
    private:
        mat k = make_xd_mat(6,6);
        mat N = make_xd_mat(2,6);
        mat B = make_xd_mat(2,6);
        std::vector<int> dof_map = {0, 2, 5};
    public:
        mat const get_k() const {return k;}
        mat const get_N() const {return N;}
        mat const get_B() const {return B;}
        std::vector<int> const get_dof_map() const {return dof_map;};
        void calc_N(real x, real L);
        void calc_B(real x, real L);
        void calc_k(real L, BasicSection& sec);
        /**
         * @brief calculates the element material stiffness in-place for the element.
         * @details for basic beam-column element, the material stiffness is a simple predefine matrix that uses length, EA, and EI.
         * This function takes the element material stiffness by reference and places the relevant components there.
         * 
         * @param L beam length.
         * @param sec section object containing methods to retrieve axial and bending rigidities EA and EI.
         * @param k local element material stiffness matrix.
         */
        void calc_elem_mat_stiffness(real& L, BasicSection& sec, mat& k);
        
        /**
         * @brief calculates the element geometric stiffness in-place for the element.
         * @details for basic beam-column element, the geometric stiffness used here is that related to the
         * axial force of the Hermitian beam-column element. Referenced from Chapter 16 of Felippa's Nonlinear FEA notes.
         * Equation 16.26.
         * 
         * @param L beam length
         * @param P axial force in the beam-column.
         * @param k_g local element geometric stiffness matrix.
         */
        void calc_elem_geom_stiffness(real& L, real& P, mat& k_g);
};

#endif