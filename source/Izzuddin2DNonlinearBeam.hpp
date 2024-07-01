/**
 * @file Izzuddin2DNonlinearBeam.hpp
 * @brief The nonlinear 2D beam-column element introduced by Bassam Izzuddin in his NLSA notes.
 * 
 */

#ifndef IZZUDDIN_2D_NONLINEAR_BEAM_HPP
#define IZZUDDIN_2D_NONLINEAR_BEAM_HPP
#include "beam_element.hpp"
#include "NonlinearTransform.hpp"
#include "main.hpp"

class Izzuddin2DNonlinearBeam : public Basic2DBeamElement
{
    protected:
        // vec deformational_displacements = make_xd_vec(3); /**< local deformational-displacements for all freedoms which are \f$ [\delta \theta_1 \theta_2]^T \f$.*/
        NonlinearTransform corot_transform;
        
        
        
        
        // std::vector<spnz> global_R_triplets; /**< triplet vector for global resistance forces \f$\boldsymbol{R}\f$.*/
        
        mat local_mat_stiffness = make_xd_mat(3,3); /**< local element 3x3 material stiffness matrix.*/
        mat local_geom_stiffness = make_xd_mat(3,3); /**< local element 3x3 geometric stiffness matrix.*/
        mat local_tangent_stiffness = make_xd_mat(3,3); /**< local element 3x3 tangent stiffness matrix.*/
        mat external_geom_stiffness = make_xd_mat(12,12); /**< Geometric stiffness matrix contribution to global geometric stiffness due to \f$ \partial ^2\boldsymbol{d} /\partial \boldsymbol{U}^2\f$.*/
        vec local_d = make_xd_vec(3); /**< local deformational displacements \f$\boldsymbol{d} = [\Delta,    \theta_1,   \theta_2]^T\f$.*/
        vec local_f = make_xd_vec(3); /**< local nodal-forces corresponding to deformational displacements \f$ \boldsymbol{f} = [F, M_1, M_2]^T\f$.*/
    public:
        Izzuddin2DNonlinearBeam();
        
        template<typename Container>
        Izzuddin2DNonlinearBeam(int given_id, Container& in_nodes) {
            initialise(given_id, in_nodes);
        }

        template<typename Container>
        void initialise(int given_id, Container& in_nodes) {
            if (std::size(in_nodes) != nnodes)
            {
                std::cout << "Incorrect number of nodes passed to create element " << id << std::endl;
                std::cout << "Received " << std::size(in_nodes) << " but expected " << nnodes << std::endl; 
                std::exit(1);
            }
            id = given_id;
            nodes[0] = in_nodes[0];
            nodes[1] = in_nodes[1];
            for (auto node : in_nodes) {  
                node->add_connected_element(id);
            }
            corot_transform.initialise(nodes);
            calc_local_constitutive_mat();
            update_state();
        }

        /**
         * @brief calculates the material stiffness matrix using the shape-function function \ref ShapeFunction::calc_elem_mat_stiffness.
         * 
         */
        void calc_mat_stiffness()
        {
            real EI = section.get_E()*section.get_I();
            real EA = section.get_E()*section.get_A();
            // real F = local_f(0);
            // real delta = local_d(0);
            real theta1 = local_d(1);
            real theta2 = local_d(2);
            real L0 = corot_transform.get_L0();
            vec V = make_xd_vec(3);
            V(0) = 1/L0;
            V(1) = 2*theta1/15 - theta2/30;
            V(2) = -theta1/30 + 2*theta2/15;
            local_mat_stiffness.setZero();
    
            local_mat_stiffness(1,1) = 4*EI/L0;
            local_mat_stiffness(2,2) = 4*EI/L0;
            local_mat_stiffness(1,2) = 2*EI/L0;
            local_mat_stiffness(2,1) = 2*EI/L0;
            local_mat_stiffness = local_mat_stiffness + EA*L0*V*V.transpose();


        }
        /**
         * @brief calculates the geometric stiffness matrix using the shape-function function \ref ShapeFunction::calc_elem_geom_stiffness.
         * 
         */
        void calc_geom_stiffness()
        {
            real F = local_f(0);
            real L0 = corot_transform.get_L0();
            local_geom_stiffness.setZero();
            // this is for geometric stiffness!!
            local_geom_stiffness(1,1) = 4*F*L0/30;
            local_geom_stiffness(2,2) = 4*F*L0/30;
            local_geom_stiffness(1,2) = -F*L0/30;
            local_geom_stiffness(2,1) = -F*L0/30;
        }

        /**
         * @brief Creates the global contributions by transforming the tangent matrix and also adds the external stiffness contributions from Felippa (16.13), (16.14), and (16.15) and Izzuddin (10.c) and (10.d).
         */
        void calc_elem_global_stiffness()
        {
            elem_global_stiffness = corot_transform.get_nl_T().transpose()*local_tangent_stiffness*corot_transform.get_nl_T() + external_geom_stiffness;
        }

        /**
         * @brief updates the deformational displacements using the nonlinear transformation object;
         * 
         */
        void calc_d_from_U()
        {
            corot_transform.calc_deformational_displacements(local_d);
        }
        
        /**
         * @brief calculates \ref local_f based on (6.b) from Izzuddin.
         * 
         */
        void calc_local_f()
        {
            real EI = section.get_E()*section.get_I();
            real EA = section.get_E()*section.get_A();
            real F = local_f(0);
            real delta = local_d(0);
            real theta1 = local_d(1);
            real theta2 = local_d(2);
            real L0 = corot_transform.get_L0();
            local_f(0) = EA*((delta/L0) + (2*theta1*theta1 - theta1*theta2 + 2*theta2*theta2)/30);
            local_f(1) = ((4*EI/L0) + (2*F*L0/15))*theta1 + ((2*EI/L0) - F*L0/30)*theta2;
            local_f(2) = ((2*EI/L0) - F*L0/30)*theta1 + ((4*EI/L0) + (2*F*L0/15))*theta2;
        }
        /**
         * @brief calculates strains based on (4.b) and (4.c) from Izzuddin. x is taken at midpoint of element.
         * 
         */
        void calc_eps()
        {
            real delta = local_d(0);
            real theta1 = local_d(1);
            real theta2 = local_d(2);
            real L0 = corot_transform.get_L0();
            real x = 0.5*L0;
            local_eps(0) = ((delta/L0) + (2*theta1*theta1 - theta1*theta2 + 2*theta2*theta2)/30);
            local_eps(1) = ((-4/L0) + 6*x/(L0*L0))*theta1 + ((-2/L0) + 6*x/(L0*L0))*theta2;
        }
        
        // void calc_stresses();<<< same as the regular linear beam. D*eps.
        /**
         * @brief calculates the global resistance forces contributions of this element from \f$\boldsymbol{R}^e = \frac{\partial \boldsymbol{d}}{\partial \boldsymbol{U}}\boldsymbol{f}\f$
         * 
         */
        void calc_element_global_resistance_forces()
        {
            element_resistance_forces = corot_transform.get_nl_T().transpose()*local_f;
        }
        /**
         * @brief calculates the direct external geometric stiffness contributions of the element.
         * @details
         * equation (10.c) \f$ \frac{\partial ^2\theta}{\partial \boldsymbol{U}^2} =  * \begin{bmatrix} NA & \bf{0} & \bf{1} & \bf{2} & \bf{3} & \bf{4} & \bf{5} & \bf{6} & \bf{7} & \bf{8} & \bf{9} & \bf{10} & \bf{11} \\ \bf{0} & -g_1 & 0 & g_2 & 0 & 0 & 0 & g_1 & 0 & -g_2 & 0 & 0 & 0 \\ \bf{1} & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\ \bf{2} & g_2 & 0 & g_1 & 0 & 0 & 0 & -g_2 & 0 & -g_1 & 0 & 0 & 0  \\ \bf{3} & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\ \bf{4} & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\ \bf{5} & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\ \bf{6} & g_1 & 0 & -g_2 & 0 & 0 & 0 & -g_1 & 0 & g_2 & 0 & 0 & 0 \\ \bf{7} & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\ \bf{8} & -g_2 & 0 & -g_1 & 0 & 0 & 0 & g_2 & 0 & g_1 & 0 & 0 & 0  \\ \bf{9} & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\ \bf{10} & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\ \bf{11} & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \end{bmatrix}\f$
         * 
         * equation (10.d) \f$ \frac{\partial ^2\Delta}{\partial \boldsymbol{U}^2} =  * \begin{bmatrix} NA & \bf{0} & \bf{1} & \bf{2} & \bf{3} & \bf{4} & \bf{5} & \bf{6} & \bf{7} & \bf{8} & \bf{9} & \bf{10} & \bf{11} \\ \bf{0} & g_5 & 0 & -g_4 & 0 & 0 & 0 & -g_5 & 0 & g_4 & 0 & 0 & 0 \\ \bf{1} & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\ \bf{2} & -g_4 & 0 & g_3 & 0 & 0 & 0 & g_4 & 0 & -g_3 & 0 & 0 & 0  \\ \bf{3} & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\ \bf{4} & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\ \bf{5} & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\ \bf{6} & -g_5 & 0 & g_4 & 0 & 0 & 0 & g_5 & 0 & -g_4 & 0 & 0 & 0 \\ \bf{7} & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\ \bf{8} & g_4 & 0 & -g_3 & 0 & 0 & 0 & -g_4 & 0 & g_3 & 0 & 0 & 0  \\ \bf{9} & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\ \bf{10} & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\ \bf{11} & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \end{bmatrix}\f$
         */
        void calc_external_geom_stiffness()
        {
            mat d2delta_du2 = make_xd_mat(12,12);
            mat d2theta_du2 = make_xd_mat(12,12);
            real g1, g2, g3, g4, g5;
            g1 = corot_transform.get_g1();
            g2 = corot_transform.get_g2();
            g3 = corot_transform.get_g3();
            g4 = corot_transform.get_g4();
            g5 = corot_transform.get_g5();

            d2delta_du2(0,0) = g5;
            d2delta_du2(0,2) = -g4;
            d2delta_du2(0,6) = -g5;
            d2delta_du2(0,8) = g4;

            d2delta_du2(2,0) = -g4;
            d2delta_du2(2,2) = g3;
            d2delta_du2(2,6) = g4;
            d2delta_du2(2,8) = -g3;

            d2delta_du2(6,0) = -g5;
            d2delta_du2(6,2) = g4;
            d2delta_du2(6,6) = g5;
            d2delta_du2(6,8) = -g4;

            d2delta_du2(8,0) = g4;
            d2delta_du2(8,2) = -g3;
            d2delta_du2(8,6) = -g4;
            d2delta_du2(8,8) = g3;

            d2theta_du2(0,0) = -g1;
            d2theta_du2(0,2) = g2;
            d2theta_du2(0,6) = g1;
            d2theta_du2(0,8) = -g2;
            
            d2theta_du2(2,0) = g2;
            d2theta_du2(2,2) = g1;
            d2theta_du2(2,6) = -g2;
            d2theta_du2(2,8) = -g1;

            d2theta_du2(6,0) = g1;
            d2theta_du2(6,2) = -g2;
            d2theta_du2(6,6) = -g1;
            d2theta_du2(6,8) = g2;

            d2theta_du2(8,0) = -g2;
            d2theta_du2(8,2) = -g1;
            d2theta_du2(8,6) = g2;
            d2theta_du2(8,8) = g1;

            external_geom_stiffness = d2delta_du2*local_f(0) + d2theta_du2*local_f(1) + d2theta_du2*local_f(2);

        }
        void calc_tangent_stiffness() {
            if (VERBOSE_STIFFNESSES)
            {
            std::cout << ">>>>>>>>>>>Izzuddin's calc_tangent_stifness<<<<<<<<<<" << std::endl;
            std::cout << "calc_tangent_stiffness::elem " << id << " mat_stiffness is " << std::endl << local_mat_stiffness << std::endl;
            std::cout << "calc_tangent_stiffness::elem " << id << " geom_stiffness is " << std::endl << local_geom_stiffness << std::endl;
            std::cout << "calc_tangent_stiffness::elem " << id << " mat_stiffness + geom_stiffness is " << std::endl << local_mat_stiffness + local_geom_stiffness << std::endl;
            }
            this->local_tangent_stiffness = this->local_mat_stiffness + this->local_geom_stiffness;
            if (VERBOSE_STIFFNESSES)
            {
            std::cout << "calc_tangent_stiffness::elem " << id << " tangent_stiffness is " << std::endl << local_tangent_stiffness << std::endl;
            }
        }

        void update_state()
        {
            if (VERBOSE_NLB)
            {
                std::cout << "element " << id << " global-level disp prior to get_U_from_nodes = " << std::endl << global_ele_U << std::endl;
                corot_transform.print_state();
            }
            
            get_U_from_nodes();
            corot_transform.update_state(global_ele_U);
            if (VERBOSE_NLB)
            {
                std::cout << "element " << id << " global-level disp after get_U_from_nodes and updating state = " << std::endl << global_ele_U << std::endl;
                corot_transform.print_state();
            }
            calc_d_from_U();
            if (VERBOSE_NLB)
            {
                std::cout << "d = " << std::endl << local_d << std::endl;
            }
            calc_local_f();
            calc_eps();
            calc_stresses();
            if (VERBOSE_NLB)
            {
                std::cout << "f = " << std::endl << local_f << std::endl;
                std::cout << "eps = " << std::endl << local_eps << std::endl;
                std::cout << "stresses = " << std::endl << local_stresses << std::endl;
            }
            calc_mat_stiffness();
            if (VERBOSE_STIFFNESSES)
                std::cout << std::endl << "--------------------------------------------------------------" << std::endl;
            std::cout << "elem " << id << " mat_stiffness is " << std::endl << local_mat_stiffness << std::endl;
            calc_geom_stiffness();
            if (VERBOSE_STIFFNESSES)
            {
            std::cout << "elem " << id << " geom_stiffness is " << std::endl << local_geom_stiffness << std::endl;
            std::cout << "elem " << id << " mat_stiffness + geom_stiffness is " << std::endl << local_mat_stiffness + local_geom_stiffness << std::endl;
            }
            calc_tangent_stiffness();
            if (VERBOSE_STIFFNESSES)
                std::cout << "elem " << id << " tan_stiffness is " << std::endl << local_tangent_stiffness << std::endl;
            calc_external_geom_stiffness();
            if (VERBOSE_STIFFNESSES)
                std::cout << "elem " << id << " external_geom_stiffness is " << std::endl << external_geom_stiffness << std::endl;
            calc_elem_global_stiffness();
            if (VERBOSE_STIFFNESSES)
                std::cout << "elem " << id << " global_stiffness contribution is " << std::endl << elem_global_stiffness << std::endl;
            populate_resistance_force_triplets();
        }
        
    
};
#endif