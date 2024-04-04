/**
 * @file NonlinearTransform.hpp
 * @brief Nonlinear transformation matrix for geometrically nonlinear analysis.
 * 
 */

#ifndef NONLINEAR_TRANSFORM_HPP
#define NONLINEAR_TRANSFORM_HPP
#include "ElementConfiguration.hpp"
/**
 * @brief NonlinearTransform object is responsible for creating the geometerically-nonlinear transformation matrix for beam-column analysis.
 * 
 * @details based on Felippa's Nonlinear Finite Element Method
 * @todo should inherit many of the functionality of the Orientation object which only needs to host an axis system and do some simple trig operations.
 */
class NonlinearTransform
{
private:
    ElementConfiguration base_configuration;
    ElementConfiguration corotated_configuration;
    ElementConfiguration current_configuration;

    real x21, y21; /**<  Geometric quantities from  from Felippa \f$ x_{21}\f$ and Izzuddin \f$ \hat{X}_E \f$ \f$ x_{21} = X_{21} + u_{x2} - u_{x1} = \hat{X}_E\f$*/
    vec* global_ele_U; /**< pointer to global nodal displacements for all freedoms of the nodes corresponding to the element. Passed from the element.*/
    real ux1, ux2, uy1, uy2, theta1, theta2; /**< Global element displacements.*/

    real L0, L; /**< L0 and L are the original length and the current length.*/
    
    real Psi; /**< angle \f$\Psi\f$ between base configuration and current (and equally corotated) configuration. \bf{Rigid body rotation!}*/
    real varphi; /**< angle \f$\varphi\f$ between global axes and base configuration.*/
    real phi; /**< summation of \ref Psi and \ref varphi - total angle between global axes and current/corotated configuration.*/
    real cos_phi, sin_phi;
    real cos_varphi, sin_varphi;
    real cos_Psi, sin_Psi;

    mat nl_T = make_xd_mat(3,12); /**< The 3x12 T transformation matrix \f$ \partial \boldsymbol{d}/\partial \boldsymbol{U}\f$.*/
public:
    void calc_global_elem_disps()
    {
        ux1 = (*global_ele_U)(0);
        ux2 = (*global_ele_U)(6);
        uy1 = (*global_ele_U)(1);
        uy2 = (*global_ele_U)(7);
        // zz-rotations.
        theta1 = (*global_ele_U)(5);
        theta2 = (*global_ele_U)(11);
    }
    /**
     * @brief calculates geometric quantities from  from Felippa \f$ x_{21}, y_{21}\f$ and Izzuddin \f$ \hat{X}_E, \hat{Y}_E\f$.
     * @details equations 8.c and 8.d from Izzuddin, equation 11.30 from Felippa.
     * 
     */
    void calc_distance_between_nodes()
    {
        x21 = base_configuration.X21 + ux2 - ux1;
        y21 = base_configuration.Y21 + uy2 - uy1;
    }
    /**
     * @brief calculates current L from equation 11.30 from Felippa: \f$ L=\sqrt{x_{21}^2 + y_{21}^2}.\f$
     * 
     */
    void calculate_L()
    {
        L = std::sqrt(x21*x21 + y21*y21);
    }
    /**
     * @brief grab L0 from the base configuration.
     * 
     */
    void initiate_L0()
    {
        L0 = base_configuration.L;
    }
    /**
     * @brief calculates the trigonometric identities based on 11.28, 11.29, 11.30, and 11.31 from Felippa.
     * 
     */
    void calc_trigonometric_identities()
    {
        cos_phi = x21/L;
        sin_phi = y21/L;
        phi = std::acos(cos_phi);

        cos_varphi = base_configuration.X21/L0;
        sin_varphi = base_configuration.Y21/L0;
        varphi = std::acos(cos_varphi);

        cos_Psi = (base_configuration.X21*x21 + base_configuration.Y21*y21)/(L*L0);
        sin_Psi = (base_configuration.X21*y21 - base_configuration.Y21*x21)/(L*L0);
        Psi = std::acos(cos_Psi);
    }
    /**
     * @brief calculations the deformational displacements based on 16.11.
     * 
     * @details we do not have vertical displacements corresponding to y because those would result in rotation
     * of the configuration which is then captured by Psi!
     * 
     * @param deformational_disp 
     */
    void calc_deformational_displacements(vec& deformational_disp)
    {
        deformational_disp(0) = L - L0;
        deformational_disp(1) = theta1 - Psi;
        deformational_disp(2) = theta2 - Psi;
    }
    /**
     * @brief  creates the T matrix which is \f$ \partial \boldsymbol{d}/\partial \boldsymbol{U}\f$.
     * 
     * @details this matrix will be 3x12. Follows 16.12 from Felippa and 9.c from Izzuddin.
     */
    void calc_T() {
        real c = std::cos(alpha);
        real s = std::sin(alpha);
        T(0,0) = c;
        T(0,2) = s;
        T(0,5) = offset;
        T(1,0) = -s;
        T(1,2) = c;
        T(2,5) = 1;
        T(3,6) = c;
        T(3,8) = s;
        T(3,11) = offset;
        T(4,6) = -s;
        T(4,8) = c;
        T(5,11) = 1;
    }

};


#endif