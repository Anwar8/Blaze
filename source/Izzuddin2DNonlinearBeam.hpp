/**
 * @file Izzuddin2DNonlinearBeam.hpp
 * @brief The nonlinear 2D beam-column element introduced by Bassam Izzuddin in his NLSA notes.
 * 
 */

#ifndef IZZUDDIN_2D_NONLINEAR_BEAM_HPP
#define IZZUDDIN_2D_NONLINEAR_BEAM_HPP
#include "beam_element.hpp"
#include "NonlinearTransform.hpp"

class Izzuddin2DNonlinearBeam : protected Basic2DBeamElement
{
    protected:
        vec deformational_displacements = make_xd_vec(3); /**< local deformational-displacements for all freedoms which are \f$ [\delta \theta_1 \theta_2]^T \f$.*/
        NonlinearTransform corot_transform;

        // mat local_mat_stiffness = make_xd_mat(3,3); /**< local element material stiffness matrix.*/
        // mat local_geom_stiffness = make_xd_mat(3,3); /**< local element geometric stiffness matrix.*/
        // mat local_tangent_stiffness = make_xd_mat(3,3); /**< local element tangent stiffness matrix.*/
        // vec local_d = make_xd_vec(3); /**< local nodal-displacements for all freedoms.*/
        vec local_f = make_xd_vec(3); /**< local nodal-forces corresponding to all freedoms.*/
    public:
        Izzuddin2DNonlinearBeam();

        // /**
        //  * @brief calculates the material stiffness matrix using the shape-function function \ref ShapeFunction::calc_elem_mat_stiffness.
        //  * 
        //  */
        // void calc_mat_stiffness();
        // /**
        //  * @brief calculates the geometric stiffness matrix using the shape-function function \ref ShapeFunction::calc_elem_geom_stiffness.
        //  * 
        //  */
        // void calc_geom_stiffness();
        // /**
        //  * @brief sums the \ref local_mat_stiffness and \ref local_geom_stiffness to create the tangent stiffness.
        //  * 
        //  */
        // void calc_tangent_stiffness() {local_tangent_stiffness = local_mat_stiffness + local_geom_stiffness;}

        /**
         * @brief Creates the global contributions by transforming the tangent matrix and also adds the external stiffness contributions from Felippa (16.13), (16.14), and (16.15) and Izzuddin (10.c) and (10.d).
         */
        void calc_elem_global_stiffness()
        {
            elem_global_stiffness = corot_transform.get_T().transpose()*local_tangent_stiffness*corot_transform.get_T();
        }

        /**
         * @brief updates the deformational 
         * 
         */
        void calc_d_from_U()
        {
            corot_transform.calc_deformational_displacements(deformational_displacements);
        }
        

        void calc_nodal_forces();

        void calc_stresses();
        void calc_global_resistance_forces();
        void calc_external_geom_stiffness();
        void update_state();
    
};
#endif