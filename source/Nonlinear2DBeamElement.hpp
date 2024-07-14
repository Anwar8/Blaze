/**
 * @file Nonlinear2DBeamElement.hpp
 * @brief definitions for the basic Euler-Bernouli nonlinear beam-column from Izzuddin's NLSA notes.
 */

#ifndef NONLINEAR_2D_BEAM_ELEMENT_HPP
#define NONLINEAR_2D_BEAM_ELEMENT_HPP

#include "BeamElementCommonInterface.hpp"

/**
 * @brief a geometrically nonlinear 2D beam element with 6 total freedoms: 1 rotation and two displacements
 * 
 * @details
 * This beam-column element has fewer freedoms than required in a 3D domain, and so the
 * transformation matrix T must map the element from 6 freedoms to the required 12 in 3D
 * domains.
 * 
 * This is a geometrically nonlinear beam-column element. It is based on the Euler-Bernouli model
 * and has a cubic shape function, and is built based on Bassam Izzuddin's NLSA notes, and also
 * on Felippa's nonlinear finite element notes.
 * 
 */
class Nonlinear2DBeamElement : public BeamElementCommonInterface {
    private:
    protected:
        /**
         * @name basic_information
         * @brief the basic data about the generic beam element.
         * all basic_information is inherited from \ref BeamElementBaseClass.
         */
        //@{
            real initial_length; //**<the initial length of the element before any deformation.*/
        //@}

        /**
         * @name beam_basic_objects
         * @brief basic objects needed by the beam-column elements. Section, shape function, transformation, etc.
         * all beam_basic_objects are inherited from \ref BeamElementBaseClass
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
        Nonlinear2DBeamElement(int given_id, Container& in_nodes) {
            initialise(given_id, in_nodes);
        }

        /**
         * @brief initialises the beam column-element with an id, nodes, and initialise any parameters needed for a new element.
         * 
         * @tparam Container any type of std container that has a std::size and built-in iterators
         * @param given_id unique identifier for the element; will be passed to the nodes
         * @param in_nodes a container of shared pointers to node objects
         */
        template<typename Container>
        void initialise(int given_id, Container& in_nodes) {
            // initialise the fundamental aspects of the element
            // -----------------------------------------------------
            elem_type = "Nonlinear_2D_EulerBernouli_beam-column"; /**< string that represents the type of the element.*/
            ndofs = 3; /**< number of freedoms at each node. 3 for this element type.*/
            nnodes = 2; /**< number of nodes. 2 nodes for this element type.*/
            set_gauss_points(); /**< set the gauss points (numbers and locations) for the element.*/
            initialise_state_containers();
            // -----------------------------------------------------

            if (std::size(in_nodes) != nnodes)
            {
                std::cout << "Incorrect number of nodes passed to create element " << id << std::endl;
                std::cout << "Received " << std::size(in_nodes) << " but expected " << nnodes << std::endl; 
                std::exit(1);
            }
            id = given_id;
            nodes.push_back(in_nodes[0]);
            nodes.push_back(in_nodes[1]);
            for (auto& node : in_nodes) {
                
                node->add_connected_element(id);
            }
            // CAUTION: copied from Izzuddin2DNonlinearBeam.hpp
            transformation.initialise(nodes);
            initial_length = transformation.get_L0(); // the initial length does not change and only needs to be calculated once.
            calc_local_constitutive_mat();
            update_state();
            // calc_T();
            // calc_length();
            // calc_local_constitutive_mat();
            // calc_stiffnesses();
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
            local_eps = make_xd_vec(2); /**< local strains. \f$ \boldsymbol{\varepsilon} = \begin{bmatrix} \varepsilon_{xx} & \kappa\end{bmatrix}^T\f$*/
            local_stresses = make_xd_vec(2); /**< local stresses.\f$ \boldsymbol{\sigma} = \begin{bmatrix} N & M \end{bmatrix}^T\f$*/
            N = make_xd_mat(2,6); /**< the shape function of the element. For this 2D element, that is 2 rows and 6 columns.*/
            B = make_xd_mat(2,6); /**< the derivative of the shape function of the element. In this case 2 rows and 6 columns.*/   
            local_constitutive_mat = make_xd_mat(2,2); /**< local constitutive matrix \f$\boldsymbol{D} = \begin{bmatrix} EA & 0 \\ 0 & EI\end{bmatrix}\f$.*/
            local_mat_stiffness = make_xd_mat(3,3); /**< local element 3x3 material stiffness matrix.*/
            local_geom_stiffness = make_xd_mat(3,3); /**< local element 3x3 geometric stiffness matrix.*/
            local_tangent_stiffness = make_xd_mat(3,3); /**< local element 3x3 tangent stiffness matrix.*/
            elem_global_stiffness = make_xd_mat(12,12); /**< the global contribution of the element - as in, tangent stiffness after transform via \f$ \boldsymbol{K}_t^e = \boldsymbol{T}^T \boldsymbol{k}_t \boldsymbol{T}\f$*/
            external_geom_stiffness = make_xd_mat(12,12); /**< Geometric stiffness matrix contribution to global geometric stiffness - a 12x12 matrix.*/
            
        }

        /**
         * @brief Set the gauss points std::vector.
         * @details Sets \ref gauss_points to an appropriate size and set of values. To be called by constructor. 
         * Not actually needed for this beam-column element without material nonlinearity, but will become necessary for materially-nonlinear versions.
         * 
         */
        void set_gauss_points() {
            gauss_points = {-0.57735, 0.57735};
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
         * @brief calculates the shape function of the beam-column element for location x.
         * @details this shape function is the same as the linear element \ref Linear2DBeamElement shape function but uses the current deformed length in stead of the initial length.
         * @warning Must call \ref calc_length before this function to ensure the length calculation is correct.
         * @param x location at which shape function is evaluated.
         */
        void calc_N(real x)
        {
            N(0,0) = 1 - (x/length);
            N(0,3) = x / length;
            N(1,1) = 1 - 3*std::pow(x/length,2) + 2*std::pow(x/length,3);
            N(1,2) = x - 2*std::pow(x,2)/length + std::pow(x/length, 2)*x;
            N(1,4) = 3*std::pow(x/length, 2) - 2*std::pow(x/length, 3);
            N(1,5) = -x*(x/length) + x * std::pow(x/length,2);
        }
        /**
         * @brief call the shape function's derivative of the shape function operation to calculate at a specific point.
         * @details this shape function derivative is the same as the linear element \ref Linear2DBeamElement shape function derivative but uses the current deformed length in stead of the initial length.
         * @warning Must call \ref calc_length before this function to ensure the length calculation is correct.
         * @param x location along beam-column element at which to calculate the derivative of the shape function.
         */
        void calc_B(real x) 
        {
            B(0,0) = -1/length;
            B(0,3) = 1/length;
            B(1,1) = -6*std::pow(1/length,2) + 12*x*std::pow(1/length,3);
            B(1,2) = - 4/length + 6*x*std::pow(1/length, 2);
            B(1,4) = 6*std::pow(1/length, 2) - 12*x*std::pow(1/length, 3);
            B(1,5) = -2/length + 6 * x* std::pow(1/length,2);
        }

        /**
         * @brief calculates the transformation matrix from the \ref NonlinearTransform \ref NonlinearTransfor::get_T function.
         * @warning: doesn not use the input variables. Just pass nothing, please, as the function needs to match the generic interface defined in \ref BeamElementBaseClass.
         * @param sec_offset is the y-axis offset of the section nodes from their actual coordinates. IGNORED.
         * @param origin_x the x-axis of the global coordinate system. IGNORED.
         */
        void calc_T(real sec_offset = 0.0, coords origin_x = {1.0, 0.0, 0.0}) {
            // orient.evaluate(nodes, sec_offset, origin_x);  
            transformation.get_T(); 
        }

        /**
         * @brief calculates local constitutive matrix from section information. For this lineaer element it is simply EA and EI along diagonals 1,1 and 2,2
         * 
         */
        void calc_local_constitutive_mat() {
            // given all constitutive mat elements are zeroed we only need to calculate the non-zero diagonal members of this element.
            local_constitutive_mat(0,0) = section.get_E()*section.get_A();
            local_constitutive_mat(1,1) = section.get_E()*section.get_I();
        }

        /**
         * @brief calculates strains based on (4.b) and (4.c) from Izzuddin. x is taken at midpoint of element.
         * @details the relationship between the \ref local_d and \ref local_eps is nonlinear, and as such we 
         * cannot simply use \f$\boldsymbol{\sigma} = \boldsymbol{B}\boldsymbol{d}\f$ which was used in \ref Linear2DBeamElement.
         */
        void calc_eps() {

            real delta = local_d(0);
            real theta1 = local_d(1);
            real theta2 = local_d(2);

            real x = 0.5*initial_length;
            local_eps(0) = ((delta/initial_length) + (2*theta1*theta1 - theta1*theta2 + 2*theta2*theta2)/30);
            local_eps(1) = ((-4/initial_length) + 6*x/(initial_length*initial_length))*theta1 + ((-2/initial_length) + 6*x/(initial_length*initial_length))*theta2;

        }

        /**
         * @brief calculates the local stresses from \f$\boldsymbol{\sigma}=\boldsymbol{D}{\boldsymbol{\varepsilon}}\f$.
         * @details despite being a nonlinear element, the stress calculation remains simply as \f$ \boldsymbol{\sigma} = \boldsymbol{D} \boldsymbol{\varepsilon}\f$
         * @warning requires \ref local_constitutive_mat and \ref local_eps to be calculated before this function is called.
         */
        void calc_stresses()  {local_stresses = local_constitutive_mat*local_eps;}
        
        /**
         * @brief  calculates \ref local_f based on (6.b) from Izzuddin noting we changed the order of the DoFs.
         * @details
         * \f$ \boldsymbol{f} = \begin{bmatrix} F \\ M_1 \\ M_2 \end{bmatrix} 
         * = \begin{bmatrix} F_i + EA\left( \frac{\Delta}{L} + \frac{2\theta_1^2 - \theta_1 \theta_2 + 2 \theta_2 ^2}{30}\right) \\
         * \left( \frac{4 EI}{L_0} + \frac{2FL_0}{15}\right) \theta_1 + \left(\frac{2EI}{L_0}  - \frac{FL_0}{30}\right)\theta_2 \\ 
         * \left( \frac{2 EI}{L_0} - \frac{FL_0}{30}\right) \theta_1 + \left(\frac{4EI}{L_0}  + \frac{2FL_0}{15}\right)\theta_2 \end{bmatrix}\f$
         */
        void calc_local_f() {

            real EI = local_constitutive_mat(1,1);
            real EA = local_constitutive_mat(0,0);

            real F = local_f(0);
            real delta = local_d(0);
            real theta1 = local_d(1);
            real theta2 = local_d(2);
            
            local_f(0) = EA*((delta/initial_length) + (2*theta1*theta1 - theta1*theta2 + 2*theta2*theta2)/30);
            local_f(1) = ((4*EI/initial_length) + (2*F*initial_length/15))*theta1 + ((2*EI/initial_length) - F*initial_length/30)*theta2;
            local_f(2) = ((2*EI/initial_length) - F*initial_length/30)*theta1 + ((4*EI/initial_length) + (2*F*initial_length/15))*theta2;

        }
                
        /**
         * @brief updates element nodal displacements, strains, stresses, element resistance forces.
         * @warning this is a purely linear element so no re-evaluation of material strength or constitutive matrix.
         * @warning This calculation is NOT done for Gauss points, but only for midpoint of element (which is preferred for output postprocessing).
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
            calc_stresses();
            calc_local_f();

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
         * @brief calculates the material stiffness matrix using (7.b) from Izzuddin's notes being careful to reorder the DoFs.
         * 
         */
        void calc_mat_stiffness() {
            real EA = local_constitutive_mat(0,0);
            real EI = local_constitutive_mat(1,1);
            
            real theta1 = local_d(1);
            real theta2 = local_d(2);
            vec V = make_xd_vec(3);

            V(0) = 1/initial_length;
            V(1) = 2*theta1/15 - theta2/30;
            V(2) = -theta1/30 + 2*theta2/15;
            local_mat_stiffness.setZero();
    
            local_mat_stiffness(1,1) = 4*EI/initial_length;
            local_mat_stiffness(2,2) = 4*EI/initial_length;
            local_mat_stiffness(1,2) = 2*EI/initial_length;
            local_mat_stiffness(2,1) = 2*EI/initial_length;
            local_mat_stiffness = local_mat_stiffness + EA*initial_length*V*V.transpose();
        }

        /**
         * @brief calculates the geometric stiffness from Izzuddin's notes equation (7.d).
         * 
         */
        void calc_geom_stiffness() 
        {
            real F = local_f(0);
            local_geom_stiffness(1,1) = 4*F*initial_length/30;
            local_geom_stiffness(2,2) = 4*F*initial_length/30;
            local_geom_stiffness(1,2) = -F*initial_length/30;
            local_geom_stiffness(2,1) = -F*initial_length/30;
        }
        /**
         * @brief calculates the direct external geometric stiffness contributions of the element.
         * @details
         * equation (10.c) \f$ \frac{\partial ^2\theta}{\partial \boldsymbol{U}^2} =  * \begin{bmatrix} NA & \bf{0} & \bf{1} & \bf{2} & \bf{3} & \bf{4} & \bf{5} & \bf{6} & \bf{7} & \bf{8} & \bf{9} & \bf{10} & \bf{11} \\ \bf{0} & -g_1 & 0 & g_2 & 0 & 0 & 0 & g_1 & 0 & -g_2 & 0 & 0 & 0 \\ \bf{1} & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\ \bf{2} & g_2 & 0 & g_1 & 0 & 0 & 0 & -g_2 & 0 & -g_1 & 0 & 0 & 0  \\ \bf{3} & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\ \bf{4} & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\ \bf{5} & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\ \bf{6} & g_1 & 0 & -g_2 & 0 & 0 & 0 & -g_1 & 0 & g_2 & 0 & 0 & 0 \\ \bf{7} & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\ \bf{8} & -g_2 & 0 & -g_1 & 0 & 0 & 0 & g_2 & 0 & g_1 & 0 & 0 & 0  \\ \bf{9} & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\ \bf{10} & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\ \bf{11} & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \end{bmatrix}\f$
         * 
         * equation (10.d) \f$ \frac{\partial ^2\Delta}{\partial \boldsymbol{U}^2} =  * \begin{bmatrix} NA & \bf{0} & \bf{1} & \bf{2} & \bf{3} & \bf{4} & \bf{5} & \bf{6} & \bf{7} & \bf{8} & \bf{9} & \bf{10} & \bf{11} \\ \bf{0} & g_5 & 0 & -g_4 & 0 & 0 & 0 & -g_5 & 0 & g_4 & 0 & 0 & 0 \\ \bf{1} & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\ \bf{2} & -g_4 & 0 & g_3 & 0 & 0 & 0 & g_4 & 0 & -g_3 & 0 & 0 & 0  \\ \bf{3} & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\ \bf{4} & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\ \bf{5} & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\ \bf{6} & -g_5 & 0 & g_4 & 0 & 0 & 0 & g_5 & 0 & -g_4 & 0 & 0 & 0 \\ \bf{7} & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\ \bf{8} & g_4 & 0 & -g_3 & 0 & 0 & 0 & -g_4 & 0 & g_3 & 0 & 0 & 0  \\ \bf{9} & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\ \bf{10} & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\ \bf{11} & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \end{bmatrix}\f$
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
