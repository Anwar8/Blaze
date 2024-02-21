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
};

#endif