/**
 * @file NonlinearTransform.hpp
 * @brief Nonlinear transformation matrix for geometrically nonlinear analysis.
 * 
 */

#ifndef NONLINEAR_TRANSFORM_HPP
#define NONLINEAR_TRANSFORM_HPP
#include "ElementConfiguration.hpp"
#include "node.hpp"
#include "main.hpp"
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
    // ElementConfiguration current_configuration; // no point to have this - not used!

    real x21, y21; /**<  Geometric quantities from  from Felippa \f$ x_{21}\f$ and Izzuddin \f$ \hat{X}_E \f$ \f$ x_{21} = X_{21} + u_{x2} - u_{x1} = \hat{X}_E\f$*/
    real ux1, ux2, uy1, uy2, theta1, theta2; /**< Global element displacements.*/

    real L0, L; /**< L0 and L are the original length and the current length.*/
    
    real psi; /**< angle \f$\psi\f$ between base configuration and current (and equally corotated) configuration. \bf{Rigid body rotation!}*/
    real varphi; /**< angle \f$\varphi\f$ between global axes and base configuration.*/
    real phi; /**< summation of \ref psi and \ref varphi - total angle between global axes and current/corotated configuration.*/
    real cos_phi, sin_phi;
    real cos_varphi, sin_varphi;
    real cos_psi, sin_psi;

    mat nl_T = make_xd_mat(3,12); /**< The 3x12 T transformation matrix \f$ \partial \boldsymbol{d}/\partial \boldsymbol{U}\f$.*/
    mat T = make_xd_mat(6,12); /**< The 6x12 T transformation matrix from (16.24).*/
public:
    
    /**
     * @brief initialise the nonlinear transform object with a direct reference to the nodes to extract their initial coordinates.
     * 
     * @param nodes 
     */
    void initialise(std::array<std::shared_ptr<Node>, 2>& nodes)
    {
        coords pt1, pt2;
        pt1 = nodes[0]->get_coords();
        pt2 = nodes[1]->get_coords();
        base_configuration.update_pts(pt1, pt2);
        initialise_L0();
    }
    void update_state(vec& global_ele_U)
    {
        extract_global_elem_disps(global_ele_U);
        calc_distance_between_nodes();
        calculate_L();
        calc_trigonometric_identities();
        calc_nl_T();
        calc_T();

    }
    void print_state()
    {
        std::cout << "NLTransform. L0 = " << L0 << ", L = " << L << std::endl;
        std::cout << "phi = " << phi << ", cos_phi = " << cos_phi << ", sin_phi = " << sin_phi << std::endl;
    }
    /**
     * @brief extracts global element displacements from a pointer to the global element displacement that lives in the beam element object.
     * @details \f$\boldsymbol{U}^e =  \begin{bmatrix}0 & 1 & 2 & 3 & 4 & 5 & 6 & 7 & 8 & 9 & 10 & 11\\ U_1^i & U_{11}^i & U_2^i & U_{22}^i & U_3^i & U_{33}^i & U_{1}^{ii} & U_{11}^{ii} & U_{2}^{ii} & U_{22}^{ii} & U_{3}^{ii} & U_{33}^{ii} \end{bmatrix}^T\f$.
     */
    void extract_global_elem_disps(vec& global_ele_U)
    {
        ux1 = global_ele_U(0);
        ux2 = global_ele_U(6);
        uy1 = global_ele_U(1);
        uy2 = global_ele_U(7);
        // zz-rotations.
        theta1 = global_ele_U(5);
        theta2 = global_ele_U(11);
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
    void initialise_L0()
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

        cos_psi = (base_configuration.X21*x21 + base_configuration.Y21*y21)/(L*L0);
        sin_psi = (base_configuration.X21*y21 - base_configuration.Y21*x21)/(L*L0);
        psi = std::acos(cos_psi);
    }
    /**
     * @brief calculations the deformational displacements based on 16.11.
     * 
     * @details we do not have vertical displacements corresponding to y because those would result in rotation
     * of the configuration which is then captured by psi!
     */
    void calc_deformational_displacements(vec& deformational_displacements)
    {
        deformational_displacements(0) = L - L0;
        deformational_displacements(1) = theta1 - psi;
        deformational_displacements(2) = theta2 - psi;
    }
    /**
     * @brief  creates the T matrix which is \f$ \partial \boldsymbol{d}/\partial \boldsymbol{U}\f$.
     * 
     * @details this matrix will be 3x12. Follows 16.12 from Felippa and 9.c from Izzuddin.
     * However, this matrix also needs to consider that we have 12 DoFs associated with the element in
     * global system as compared to the 6 in Izzuddin's and Felippa's notes. Finally, we need to be wary
     * that Izzuddin places the \f$ \Delta\f$ defromational freedom last, while Felippa places it first
     * as I would also prefer.
     * \f$ \frac{\partial \boldsymbol{d}}{\partial \boldsymbol{U}} = \begin{bmatrix} NA & \bf{0} & \bf{1} & \bf{2} & \bf{3} & \bf{4} & \bf{5} & \bf{6} & \bf{7} & \bf{8} & \bf{9} & \bf{10} & \bf{11} \\ \bf{0} & -c_{\phi} & 0 & -s_{\phi} & 0 & 0 & 0 & c_{\phi} & 0 & s_{\phi} & 0 & 0 & 0 \\ \bf{1} & -s_{\phi}/L & 0 & c_{\phi}/L & 0 & 0 & 1 & s_{\phi}/L & 0 & -c_{\phi}/L & 0 & 0 & 0 \\ \bf{2} & -s_{\phi}/L & 0 & c_{\phi}/L & 0 & 0 & 0 & s_{\phi}/L & 0 & -c_{\phi}/L & 0 & 0 & 1\end{bmatrix}\f$
     */
    void calc_nl_T() {
        nl_T(0,0) = -cos_phi;
        nl_T(0,2) = -sin_phi;
        nl_T(0,6) = cos_phi;
        nl_T(0,8) = sin_phi;

        nl_T(1,0) = -sin_phi/L;
        nl_T(1,2) = cos_phi/L;
        nl_T(1,5) = 1;
        nl_T(1,6) = sin_phi/L;
        nl_T(1,8) = -cos_phi/L;
        
        nl_T(2,0) = -sin_phi/L;
        nl_T(2,2) = cos_phi/L;
        nl_T(2,6) = sin_phi/L;
        nl_T(2,8) = -cos_phi/L;
        nl_T(2,11) = 1;
    }
    mat get_nl_T() {return nl_T;}
        /**
     * @brief 
     * creates the T matrix
     * 
     * @details
     * the T matrix created is a 6x12 matrix with members only in the first
     * 6 rows.
     */
    void calc_T() {
        real c = cos_phi;
        real s = sin_phi;
        T(0,0) = c;
        T(0,2) = s;
        T(1,0) = -s;
        T(1,2) = c;
        T(2,5) = 1;
        T(3,6) = c;
        T(3,8) = s;
        T(4,6) = -s;
        T(4,8) = c;
        T(5,11) = 1;
    }



    mat get_T() {return T;}
    real get_L() {return L;}
    real get_L0() {return L0;}

    real get_g1()
    {
        return 2*cos_phi/(L*L);
    }
    real get_g2()
    {
        return (cos_phi*cos_phi - sin_phi*sin_phi)/(L*L);
    }
    real get_g3()
    {
        return (cos_phi*cos_phi)/(L*L);
    }
    real get_g4()
    {
        return (cos_phi*sin_phi)/L;
    }
    real get_g5()
    {
        return sin_phi*sin_phi/L;
    }
};


#endif