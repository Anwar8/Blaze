/**
 * @file Nonlinear2DPlasticBeamElement.hpp
 * @brief definitions for the basic Euler-Bernouli nonlinear beam-column element with spreading plasticity from Izzuddin's NLSA notes.
 */

#ifndef NONLINEAR_2D_PLASTIC_BEAM_ELEMENT_HPP
#define NONLINEAR_2D_PLASTIC_BEAM_ELEMENT_HPP

#include "BeamElementCommonInterface.hpp"
#include "BeamColumnFiberSection.hpp"

/**
 * @brief a geometrically nonlinear 2D beam element with 6 total freedoms: 1 rotation and two displacements
 * and the allowance for spreading plasticity.
 * @details
 * This beam-column element has fewer freedoms than required in a 3D domain, and so the
 * transformation matrix T must map the element from 6 freedoms to the required 12 in 3D
 * domains. It can use nonlinear sections and used for spreading plasticity.
 * 
 * This is a geometrically nonlinear beam-column element. It is based on the Euler-Bernouli model
 * and has a cubic shape function, and is built based on Bassam Izzuddin's NLSA notes, and also
 * on Felippa's nonlinear finite element notes. See "Materailly Nonlinear Finite Elements for One-Dimensional
 * Structural Systems" for more details on the plasticity.
 * 
 */
class Nonlinear2DPlasticBeamElement : public BeamElementCommonInterface<BeamColumnFiberSection> {
    private:
    protected:
        /**
         * @name basic_information
         * @brief the basic data about the generic beam element.
         */
        //@{
            real initial_length; //**<the initial length of the element before any deformation.*/
        //@}

        /**
         * @name beam_basic_objects
         * @brief basic objects needed by the beam-column elements. Section, shape function, transformation, etc.
         */
        //@{
        //@}

        /**
         * @name beam_state_containers
         * @brief the containers for the beam state such as displacement, strain, force, etc.
         * all beam_state_containers are inherited from \ref BeamElementBaseClass.
         */
        //@{
            mat external_geom_stiffness; /**< Geometric stiffness matrix contribution to global geometric stiffness due to \f$ \partial ^2\boldsymbol{d} /\partial \boldsymbol{U}^2\f$.*/
        //@}

    public:
    /**
     * @name element initialisers
     * @brief functions that deal with constructing and initialising the beam-column element.
     */
    //@{
        /**
         * @brief Construct a new 2D Beam Element object.
         * 
         * @tparam Container any type of std container that has a std::size and built-in iterators.
         * @param given_id unique identifier for the element; will be passed to the nodes.
         * @param in_nodes a container of shared pointers to node objects.
         */
        template<typename Container>
        Nonlinear2DPlasticBeamElement(int given_id, Container& in_nodes, BeamColumnFiberSection& sect) {
            initialise(given_id, in_nodes, sect);
        }

        /**
         * @brief initialises the beam column-element with an id, nodes, and initialise any parameters needed for a new element.
         * 
         * @tparam Container any type of std container that has a std::size and built-in iterators
         * @param given_id unique identifier for the element; will be passed to the nodes
         * @param in_nodes a container of shared pointers to node objects
         * @param sect the \ref BasicSection object that contains the material properties of the element.
         */
        template<typename Container>
        void initialise(int given_id, Container& in_nodes, BeamColumnFiberSection& sect) {
            // initialise the fundamental aspects of the element
            // -----------------------------------------------------
            elem_type = "Nonlinear_2D_EulerBernouli_Plastic_beam-column"; /**< string that represents the type of the element.*/
            ndofs = 3; /**< number of freedoms at each node. 3 for this element type.*/
            nnodes = 2; /**< number of nodes. 2 nodes for this element type.*/
            initialise_gauss_points(); /**< set the gauss points (numbers and locations) for the element.*/
            initialise_state_containers();

            if (sect.get_section_type() == Fibre)
            {
                for (int i = 0; i < gauss_points_x.size(); ++i)
                {
                    section.emplace_back(std::make_unique<BeamColumnFiberSection>(sect));
                }
            } else {
                std::cout << "Element of type " << elem_type << " only accepts section of type Fibre = 2, but got section of type: " << sect.get_section_type() << std::endl;
                exit(1);
            }
  
            
            // -----------------------------------------------------

            if (std::size(in_nodes) != nnodes)
            {
                std::cout << "Incorrect number of nodes passed to create element " << id << std::endl;
                std::cout << "Received " << std::size(in_nodes) << " but expected " << nnodes << std::endl; 
                std::exit(1);
            }
            id = given_id;
            for (auto& node: in_nodes)
            {
                nodes.emplace_back(node);
                node->add_connected_element(id);
            }
            // CAUTION: copied from Izzuddin2DNonlinearBeam.hpp
            transformation.initialise(nodes);
            initial_length = transformation.get_L0(); // the initial length does not change and only needs to be calculated once.
            update_gauss_points();
            calc_local_constitutive_mat();
            update_state();
        }

        /**
         * @brief initialises the beam state containers (matricces and vectors) to zeros of the right size for the element.
         * 
         */
        void initialise_state_containers() {
            global_ele_U = make_xd_vec(12); /**< global nodal displacements for all freedoms of the nodes corresponding to the element. 6 freedoms per node since everything is 3D, and we have two nodes for this element.*/
            local_d = make_xd_vec(3); /**< local deformational displacements \f$\boldsymbol{d} = [\Delta,    \theta_1,   \theta_2]^T\f$.*/
            local_f = make_xd_vec(3); /**< local nodal-forces corresponding to deformational displacements \f$ \boldsymbol{f} = [F, M_1, M_2]^T\f$.*/
            element_global_resistance_forces = make_xd_vec(12); /**< transformed resistance forces of the element from \ref local_f. 6 per node since we have 3D global nodes, and two nodes giving us 12 components.*/
            for (auto& gauss_point : gauss_points_x)
            {
                local_eps.emplace_back(make_xd_vec(2)); /**< local strains. \f$ \boldsymbol{\varepsilon} = \begin{bmatrix} \varepsilon_{xx} & \kappa\end{bmatrix}^T\f$*/
                local_stresses.emplace_back(make_xd_vec(2)); /**< local stresses.\f$ \boldsymbol{\sigma} = \begin{bmatrix} N & M \end{bmatrix}^T\f$*/
                N.emplace_back(make_xd_mat(2,3)); /**< the shape function of the element. For this 2D element, that is 2 rows and 6 columns.*/
                B.emplace_back(make_xd_mat(2,3)); /**< the derivative of the shape function of the element. In this case 2 rows and 6 columns.*/   
                local_constitutive_mat.emplace_back(make_xd_mat(2,2)); /**< local tangent constitutive matrix \f$\boldsymbol{D_t}\f$ which is retrieved from the fibre section.*/
            }
            local_mat_stiffness = make_xd_mat(3,3); /**< local element 3x3 material stiffness matrix.*/
            local_geom_stiffness = make_xd_mat(3,3); /**< local element 3x3 geometric stiffness matrix.*/
            local_tangent_stiffness = make_xd_mat(3,3); /**< local element 3x3 tangent stiffness matrix.*/
            elem_global_stiffness = make_xd_mat(12,12); /**< the global contribution of the element - as in, tangent stiffness after transform via \f$ \boldsymbol{K}_t^e = \boldsymbol{T}^T \boldsymbol{k}_t \boldsymbol{T}\f$*/
            external_geom_stiffness = make_xd_mat(12,12); /**< Geometric stiffness matrix contribution to global geometric stiffness - a 12x12 matrix.*/   
        }

        /**
         * @brief Set the gauss points std::vector.
         * @details Sets \ref gauss_points_x to an appropriate size and set of values. To be called by constructor. From Materially Nonlinear FEM - Izzuddin Eq (7):
         * 
         * \f$ x_{g,1} = \frac{L}{2} \left( 1 - \frac{\sqrt{3}}{3}\right) = 0.2113248654L\f$ and \f$ x_{g,2} = \frac{L}{2} \left( 1 + \frac{\sqrt{3}}{3}\right) = 0.78867513459 L\f$
         * 
         * \f$ w_{g,1} = L/2 \f$ and \f$ w_{g,2} = L/2\f$
         * 
         */
        void initialise_gauss_points() {
            gauss_points_x = {0.2113248654, 0.78867513459};
            gauss_points_w = {0.5, 0.5};
        }
        /**
         * @brief Updates Gauss points after the length of the element is known.
         * @details Multiplies \ref gauss_points_x and \ref gauss_points_w by the length of the element which is now known.
         * 
         */
        virtual void update_gauss_points()
        {
            for (auto& pt_x : gauss_points_x)
            {
                pt_x *= initial_length;
            }
            for (auto& pt_w : gauss_points_w)
            {
                pt_w *= initial_length;
            }
        }
    //@}

    /**
     * @name element property calculation functions
     * @brief functions that are used to calculate element properties such as stress, shape, length, and strain.
     */
    //@{
        /**
         * @brief retrieves the length of the element using the transformation object which is already initialised and calculates the current actual deformed length rather than the initial length.
         * 
         */
        void calc_length() {
            length = transformation.get_L();
        }
        /**
         * @brief The shape functions for the element are not needed for this element, and are not used. I might want to revisit this later.
         * @details The displacement and rotation fields can be interpolated from Izzuddin's notes equation (2), and its derivative:
         * \f$ v(x) = \left[x - \frac{2 x^2}{L_0} + \frac{x^3}{L_0 ^2} \right] \theta_1 + \left[- \frac{x^2}{L_0} + \frac{x^3}{L_0 ^2} \right] \theta_2\f$ 
         * \f$ \frac{d v(x)}{dx} = \left[1 - \frac{4 x}{L_0} + \frac{3 x^2}{L_0 ^2} \right] \theta_1 + \left[- \frac{2 x}{L_0} + \frac{3 x^2}{L_0 ^2} \right] \theta_2\f$ 
         * I assume I can easily find \f$ u(x)\f$ from my self-assumed \f$ u(x) = \frac{\Delta}{L_0} x \f$
         */
        void calc_N( )
        {
            // do nothing.
            std::cout << "WARNING: this function from Nonlinear2DBeamElement does not do anything at the moment." << std::endl;
        }
        /**
         * @brief call the shape function's derivative of the shape function operation to calculate at the Gauss points.
         * @details This shape function derivative corresponds to Izzudin's notes equation (5.b)
         * \f$\boldsymbol{B} = \frac{\partial \boldsymbol{\varepsilon}}{\partial \boldsymbol{d}^T} = \begin{bmatrix} \frac{1}{L_0} & \frac{2\theta_1}{15} - \frac{\theta_2}{30} & -\frac{\theta_1}{30} + \frac{2\theta_2}{15} \\
         * 0 & -\frac{4}{L_0} + \frac{6x}{L_0 ^2} & -\frac{2}{L_0} + \frac{6x}{L_0 ^2}\end{bmatrix} \f$ 
         * @param x location along beam-column element at which to calculate the derivative of the shape function.
         */
        void calc_B( ) 
        {
            real theta_1 = local_d(1);
            real theta_2 = local_d(2);
            for (int i = 0; i < gauss_points_x.size(); ++i)
            {
                B[i](0,0) = -1/initial_length;
                B[i](1,0) = 0;

                B[i](0,1) = 2*theta_1/15 - theta_2/30;
                B[i](1,1) = -4/initial_length + 6*gauss_points_x[i]/(initial_length*initial_length);

                B[i](0,2) = -theta_1/30 + 2*theta_2/15;
                B[i](1,2) = -2/initial_length + 6*gauss_points_x[i]/(initial_length*initial_length);
            } 
        }

        /**
         * @brief calculates the transformation matrix from the \ref NonlinearTransform \ref NonlinearTransfor::get_T function.
         * @warning doesn not use the input variables. Just pass nothing, please, as the function needs to match the generic interface defined in \ref BeamElementBaseClass.
         * @warning This function is not used in this nonlinear element.
         * @param sec_offset is the y-axis offset of the section nodes from their actual coordinates. IGNORED.
         * @param origin_x the x-axis of the global coordinate system. IGNORED.
         */
        void calc_T(real sec_offset = 0.0, coords origin_x = {1.0, 0.0, 0.0}) {
            // orient.evaluate(nodes, sec_offset, origin_x);  
            transformation.calc_T(); 
        }

        /**
         * @brief calculates local constitutive matrix from section information - retrieves \f$ \boldsymbol{D}_t \f$ from section after updating section state.
         * 
         */
        void calc_local_constitutive_mat() {
            // given all constitutive mat elements are zeroed we only need to calculate the non-zero diagonal members of this element.
            for (int i = 0; i < section.size(); ++i)
            {
                section[i]->update_section_state(local_eps[i]);
            }

            for (int i = 0; i < gauss_points_x.size(); ++i)
            {
                local_constitutive_mat[i] = section[i]->get_D_t();
            }
        }

        /**
         * @brief calculates strains based on (4.b) and (4.c) from Izzuddin. This is done per Gauss point, which in this case is just at midpoint of element.
         * @details the relationship between the \ref local_d and \ref local_eps is nonlinear, and as such we 
         * cannot simply use \f$\boldsymbol{\sigma} = \boldsymbol{B}\boldsymbol{d}\f$ which was used in \ref Linear2DBeamElement.
         */
        void calc_eps() {

            real delta = local_d(0);
            real theta1 = local_d(1);
            real theta2 = local_d(2);

            for (int i = 0; i < gauss_points_x.size(); ++i)
            {
                local_eps[i](0) = (delta/initial_length) + (2*theta1*theta1 - theta1*theta2 + 2*theta2*theta2)/30;
                local_eps[i](1) = ((-4/initial_length) + 6*gauss_points_x[i]/(initial_length*initial_length))*theta1 + ((-2/initial_length) + 6*gauss_points_x[i]/(initial_length*initial_length))*theta2;
            }
        }

        /**
         * @brief calculates the local stresses from \f$\boldsymbol{\sigma}=\boldsymbol{D}_t{\boldsymbol{\varepsilon}}\f$.
         * @details despite being a nonlinear element, the stress calculation remains simply as \f$ \boldsymbol{\sigma} = \boldsymbol{D}_t \boldsymbol{\varepsilon}\f$
         */
        void calc_stresses()  {
            for (int i = 0; i < gauss_points_x.size(); ++i)
            {
                local_stresses[i] = local_constitutive_mat[i]*local_eps[i];
            }
        }
        
        /**
         * @brief  calculates \ref local_f based on (6.b) from Izzuddin noting we changed the order of the DoFs, and that we are including nonlinearity.
         * @details
         * \f$ \boldsymbol{f} = \begin{bmatrix} F \\ M_1 \\ M_2 \end{bmatrix} 
         * = \sum_{i=1}^n w_{g,i} \left(\boldsymbol{B}^T \boldsymbol{\sigma}\right)|_{x=x_{g,i}}
         * \f$
         */
        void calc_local_f() 
        {
            real F = local_f(0);
            real delta = local_d(0);
            real theta1 = local_d(1);
            real theta2 = local_d(2);

            local_f.setZero();

            for (int i = 0; i < gauss_points_w.size(); ++i)
            {
                local_f += gauss_points_w[i] * B[i].transpose() * local_stresses[i];
            }
        }
                
        /**
         * @brief updates element nodal displacements, strains, constitutive relationship, stresses, element resistance forces.
         */
        void update_state() 
        {
            get_U_from_nodes();
            // the corotational transformation actually really cares about the global displacements, so we need to update the global displacements first.
            transformation.update_state(global_ele_U);
            // need to retrieve local displacement from the global displacements, which in our case actually uses the transformation object as the relationship is nonlinear!
            calc_d_from_U();

            // calculating element strain, stress, and resistance forces local_f depends on local displacement d in a nonlinear way.
            calc_eps();
            calc_local_constitutive_mat(); // need to update the constitutive matrix each iteration due to material nonlinearity.
            calc_stresses();
            calc_local_f(); // numerical integration required to calculate local_f.

            // stiffness calculation must come AFTER stress calculation as stiffness may depend on stress.
            calc_stiffnesses();
            
            // these local element nodal forces are transformed to global nodal forces.
            calc_element_global_resistance_forces();
            // finally, we map back the nodal forces to the force triplets that will be used to populate the global force/resistance vector.
            populate_resistance_force_triplets();
        } 
    //@}

    /**
     * @name stiffness matrix functions
     * @brief functions that deal with generating and evaluating the different stiffness matrices.
     */
    //@{
        /**
         * @brief calculates the material stiffness matrix using (5.e) from Izzuddin's notes being careful to reorder the DoFs.
         * @details 
         * \f$ \boldsymbol{k}_e = \sum_{i=1}^n w_{g,i} \left(\boldsymbol{B}^T \boldsymbol{D}_t \boldsymbol{B} \right)|_{x=x_{g,i}} \f$ 
         */
        void calc_mat_stiffness() {
            local_mat_stiffness.setZero();

            for (int i = 0; i < gauss_points_w.size(); ++i)
            {
                local_mat_stiffness += gauss_points_w[i]* (B[i].transpose() * local_constitutive_mat[i] * B[i]);
            }
        }

        /**
         * @brief calculates the geometric stiffness from Izzuddin's notes equation (5.e) where we note the absence of L due to numerical integration.
         * @details this part of the calculation is actually identical to before and does not need to be changed, but it was anyway to reflect selection
         * of Gauss points and be consistent with the theory. That is:
         * 
         * \f$ w_{g,i} = L/2\f$ 
         * 
         * and so:
         * 
         * \f$ \boldsymbol{k}_g = \frac{\partial ^2 \varepsilon_c}{\partial \boldsymbol{d} \partial \boldsymbol{d}^T}= 2 * w_{g,i} * \frac{1}{30} \begin{bmatrix}0 & 0 & 0 \\ 0 & 4 & -1 \\ 0 & -1 & 4 \end{bmatrix}\f$
         * 
         */
        void calc_geom_stiffness() 
        {
            local_geom_stiffness.setZero();
            real F = local_f(0);
            for (int i = 0; i < gauss_points_w.size(); ++i)
            {
                local_geom_stiffness(1,1) += gauss_points_w[i]*4*F*initial_length/30;
                local_geom_stiffness(2,2) += gauss_points_w[i]*4*F*initial_length/30;
                local_geom_stiffness(1,2) += gauss_points_w[i]*-F*initial_length/30;
                local_geom_stiffness(2,1) += gauss_points_w[i]*-F*initial_length/30;
            }
        }
        /**
         * @brief calculates the direct external geometric stiffness contributions of the element.
         * @details
         * equation (10.c) \f$ \frac{\partial ^2\theta}{\partial \boldsymbol{U} \partial \boldsymbol{U}^T} =  * \begin{bmatrix} NA & \bf{0} & \bf{1} & \bf{2} & \bf{3} & \bf{4} & \bf{5} & \bf{6} & \bf{7} & \bf{8} & \bf{9} & \bf{10} & \bf{11} \\ \bf{0} & -g_1 & 0 & g_2 & 0 & 0 & 0 & g_1 & 0 & -g_2 & 0 & 0 & 0 \\ \bf{1} & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\ \bf{2} & g_2 & 0 & g_1 & 0 & 0 & 0 & -g_2 & 0 & -g_1 & 0 & 0 & 0  \\ \bf{3} & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\ \bf{4} & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\ \bf{5} & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\ \bf{6} & g_1 & 0 & -g_2 & 0 & 0 & 0 & -g_1 & 0 & g_2 & 0 & 0 & 0 \\ \bf{7} & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\ \bf{8} & -g_2 & 0 & -g_1 & 0 & 0 & 0 & g_2 & 0 & g_1 & 0 & 0 & 0  \\ \bf{9} & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\ \bf{10} & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\ \bf{11} & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \end{bmatrix}\f$
         * 
         * equation (10.d) \f$ \frac{\partial ^2\Delta}{\partial \boldsymbol{U} \partial \boldsymbol{U}^T} =  * \begin{bmatrix} NA & \bf{0} & \bf{1} & \bf{2} & \bf{3} & \bf{4} & \bf{5} & \bf{6} & \bf{7} & \bf{8} & \bf{9} & \bf{10} & \bf{11} \\ \bf{0} & g_5 & 0 & -g_4 & 0 & 0 & 0 & -g_5 & 0 & g_4 & 0 & 0 & 0 \\ \bf{1} & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\ \bf{2} & -g_4 & 0 & g_3 & 0 & 0 & 0 & g_4 & 0 & -g_3 & 0 & 0 & 0  \\ \bf{3} & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\ \bf{4} & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\ \bf{5} & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\ \bf{6} & -g_5 & 0 & g_4 & 0 & 0 & 0 & g_5 & 0 & -g_4 & 0 & 0 & 0 \\ \bf{7} & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\ \bf{8} & g_4 & 0 & -g_3 & 0 & 0 & 0 & -g_4 & 0 & g_3 & 0 & 0 & 0  \\ \bf{9} & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\ \bf{10} & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\ \bf{11} & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \end{bmatrix}\f$
         * 
         * 
         * Note that the DoFs are originally \f$ \boldsymbol{U} = \begin{bmatrix} U_1 & V_1 & \alpha_1 & U_2 &  V_2 & \alpha_2 \end{bmatrix}^T\f$ 
         * 
         * The matrices have been reordered to correspond to: \f$ \boldsymbol{U} = \begin{bmatrix} {U_1}^i & {U_{11}}^i=0 & {U_2}^i & {U_{22}}^i = 0 & {U_3}^i = 0 & {U_{33}}^i &  {U_1}^{ii} & {U_{11}}^{ii}=0 & {U_2}^{ii} & {U_{22}}^{ii} = 0 & {U_3}^{ii} = 0 & {U_{33}}^{ii}  \end{bmatrix}^T\f$ where the superscript i indicates the node number.
         */
        void calc_external_geom_stiffness()
        {
            mat d2delta_du2 = make_xd_mat(12,12);
            mat d2theta_du2 = make_xd_mat(12,12);
            real g1, g2, g3, g4, g5;
            g1 = transformation.get_g1();
            g2 = transformation.get_g2();
            g3 = transformation.get_g3();
            g4 = transformation.get_g4();
            g5 = transformation.get_g5();
            /** \f$ \frac{\partial^2 \Delta}{\partial \boldsymbol{U} \partial \boldsymbol{U}^T}\f$*/
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
            /** \f$ \frac{\partial^2 \theta_1}{\partial \boldsymbol{U} \partial \boldsymbol{U}^T} = \frac{\partial^2 \theta_2}{\partial \boldsymbol{U} \partial \boldsymbol{U}^T}\f$*/
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
        /**
         * @brief sums the \ref local_mat_stiffness and \ref local_geom_stiffness to create the tangent stiffness.
         * 
         */
        void calc_tangent_stiffness() {
            local_tangent_stiffness = local_mat_stiffness + local_geom_stiffness;
        }

        /**
         * @brief calculates the global contributions of the tangent stiffness to the global stiffness matrix as a 6x6 matrix from \f$ \boldsymbol{T}^T \boldsymbol{k}_{tan} \boldsymbol{T} \f$.
         * 
         */
        void calc_elem_global_stiffness()
        {
            elem_global_stiffness = transformation.get_nl_T().transpose()*local_tangent_stiffness*transformation.get_nl_T() + external_geom_stiffness;
        }
    //@}

    /**
     * @name logging functions
     * @brief functions used for logging output to stream - mostly for debugging.
     * All logging functions are inherited from \ref BeamElementCommonInterface.
     */
    //@{
    //@}

    /**
     * @name local-global mapping function
     * @brief functions that deal mapping degrees of freedoms between local element level and global matrices.
     */
    //@{
        /**
         * @brief updates the deformational displacements using the nonlinear transformation object;
         * @details This is based on 16.11 from Felippa's notes; refer to \ref NonlinearTransform::calc_deformational_displacements 
         * for more details, but keep in mind the deformational displacements are: \f$ \boldsymbol{d} = \begin{bmatrix} \Delta & \theta_1 & \theta_2 \end{bmatrix}^T\f$
         * 
         */
        void calc_d_from_U() 
        {
            transformation.calc_deformational_displacements(local_d);
        }

        /**
         * @brief calculates the global resistance forces contributions of this element from \f$\boldsymbol{R}^e = \left[\frac{\partial \boldsymbol{d}}{\partial \boldsymbol{U}}\right]^T\boldsymbol{f}\f$
         */
        void calc_element_global_resistance_forces() {
            element_global_resistance_forces = transformation.get_nl_T().transpose()*local_f;
        }
    //@}
    /**
     * @name setter functions
     * @brief functions that set protected variables
     * @details all setter functions are inherited from \ref BeamElementCommonInterface.
     */
    //@{
    //@}

    /**
     * @name getter functions
     * @brief functions that retrieve protected variables
     * @details all getter functions are inherited from \ref BeamElementCommonInterface.
     */
    //@{
    //@}
};

#endif
