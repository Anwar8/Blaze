/**
 * @file Linear2DBeamElement.hpp
 * @brief definitions for the basic Euler-Bernouli linear beam-column.
 */

#ifndef LINEAR_2D_BEAM_ELEMENT_HPP
#define LINEAR_2D_BEAM_ELEMENT_HPP

// #include "BeamElementBaseClass.hpp"
#include "BeamElementCommonInterface.hpp"

/**
 * @brief a 2D beam element with 6 total freedoms: 1 rotation and two displacements
 * 
 * @details
 * This beam-column element has fewer freedoms than required in a 3D domain, and so the
 * transformation matrix T must map the element from 6 freedoms to the required 12 in 3D
 * domains.
 * 
 * This is a linear-only beam-column element with no geometric stiffness. It is based on 
 * the Euler-Bernouli model and has a cubic shape function. The material stiffness should
 * be identical to the Timoshenko beam-column element with the residual bending flexibility
 * correction as noted in Felippa's notes (11-18). 
 * 
 */
class Linear2DBeamElement : public BeamElementCommonInterface {
    private:
    protected:
        /**
         * @name basic_information
         * @brief the basic data about the generic beam element.
         * all basic_information is inherited from \ref BeamElementBaseClass.
         */
        //@{      
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
        template<typename Container, typename SectionType>
        Linear2DBeamElement(int given_id, Container& in_nodes, SectionType& sect) {
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
        template<typename Container, typename SectionType>
        void initialise(int given_id, Container& in_nodes, SectionType& sect) {
            // initialise the fundamental aspects of the element
            // -----------------------------------------------------
            elem_type = "2D_EulerBernouli_beam-column"; /**< string that represents the type of the element.*/
            ndofs = 3; /**< number of freedoms at each node. 3 for this element type.*/
            nnodes = 2; /**< number of nodes. 2 nodes for this element type.*/
            initialise_gauss_points(); /**< set the gauss points (numbers and locations) for the element.*/
            initialise_state_containers();
            if (sect.get_section_type() == Basic)
            {
                section.emplace_back(std::make_unique<SectionType>(sect));
            } else {
                std::cout << "Element of type " << elem_type << " only accepts section of type Basic = 1, but got section of type: " << sect.get_section_type() << std::endl;
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
            nodes.push_back(in_nodes[0]);
            nodes.push_back(in_nodes[1]);
            for (auto& node : in_nodes) {
                
                node->add_connected_element(id);
            }
            calc_T();
            calc_length();
            update_gauss_points();
            calc_local_constitutive_mat();
            calc_stiffnesses();
        }

        /**
         * @brief initialises the beam state containers (matricces and vectors) to zeros of the right size for the element.
         * 
         */
        void initialise_state_containers() {
            global_ele_U = make_xd_vec(12); /**< global nodal displacements for all freedoms of the nodes corresponding to the element. 6 freedoms per node since everything is 3D, and we have two nodes for this element.*/
            local_d = make_xd_vec(6); /**< local nodal-displacements for all freedoms. 3 DOFs per node.*/
            local_f = make_xd_vec(6); /**< local nodal-forces corresponding to all freedoms. 3 element forces per node (Fx, Fy, and Mzz).*/
            element_global_resistance_forces = make_xd_vec(12); /**< transformed resistance forces of the element from \ref local_f. 6 per node since we have 3D global nodes, and two nodes giving us 12 components.*/
            for (auto& gauss_point : gauss_points_x)
            {
                local_eps.emplace_back(make_xd_vec(2)); /**< local strains. \f$ \boldsymbol{\varepsilon} = \begin{bmatrix} \varepsilon_{xx} & \kappa\end{bmatrix}^T\f$*/
                local_stresses.emplace_back(make_xd_vec(2)); /**< local stresses.\f$ \boldsymbol{\sigma} = \begin{bmatrix} N & M \end{bmatrix}^T\f$*/
                N.emplace_back(make_xd_mat(2,6)); /**< the shape function of the element. For this 2D element, that is 2 rows and 6 columns.*/
                B.emplace_back(make_xd_mat(2,6)); /**< the derivative of the shape function of the element. In this case 2 rows and 6 columns.*/
                local_constitutive_mat.emplace_back(make_xd_mat(2,2)); /**< local constitutive matrix \f$\boldsymbol{D} = \begin{bmatrix} EA & 0 \\ 0 & EI\end{bmatrix}\f$.*/
            }
            local_mat_stiffness = make_xd_mat(6,6); /**< local element material stiffness matrix.*/
            local_geom_stiffness = make_xd_mat(6,6); /**< local element geometric stiffness matrix.*/
            local_tangent_stiffness = make_xd_mat(6,6); /**< local element tangent stiffness matrix.*/
            elem_global_stiffness = make_xd_mat(12,12); /**< the global contribution of the element - as in, tangent stiffness after transform via \f$ \boldsymbol{K}_t^e = \boldsymbol{T}^T \boldsymbol{k}_t \boldsymbol{T}\f$*/
        }

        /**
         * @brief Set the gauss points std::vector.
         * @details Sets \ref gauss_points_x to an appropriate size and set of values. To be called by constructor. 
         * Not actually needed for this beam-column element without material nonlinearity, but will become necessary for materially-nonlinear versions.
         * 
         */
        void initialise_gauss_points() {
            gauss_points_x = {0.5};
            gauss_points_w = {1.0};
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
                pt_x *= length;
            }
            for (auto& pt_w : gauss_points_w)
            {
                pt_w *= length;
            }
        }
    //@}

    /**
     * @name element property calculation functions
     * @brief functions that are used to calculate element properties such as stress, shape, length, and strain.
     */
    //@{
        /**
         * @brief retrieves the length of the element using the orientation object since this is a linear element.
         * 
         */
        void calc_length() {
            length = orient.get_length();
        }
        /**
         * @brief calculates the shape function of the beam-column element for location of Gauss points.
         */
        void calc_N()
        {            
            for (int i = 0; i < gauss_points_x.size(); ++i)
            {
                N[i](0,0) = 1 - (gauss_points_x[i]/length);
                N[i](0,3) = gauss_points_x[i] / length;
                N[i](1,1) = 1 - 3*std::pow(gauss_points_x[i]/length,2) + 2*std::pow(gauss_points_x[i]/length,3);
                N[i](1,2) = gauss_points_x[i] - 2*std::pow(gauss_points_x[i],2)/length + std::pow(gauss_points_x[i]/length, 2)*gauss_points_x[i];
                N[i](1,4) = 3*std::pow(gauss_points_x[i]/length, 2) - 2*std::pow(gauss_points_x[i]/length, 3);
                N[i](1,5) = -gauss_points_x[i]*(gauss_points_x[i]/length) + gauss_points_x[i] * std::pow(gauss_points_x[i]/length,2);
            }
        }
        /**
         * @brief call the shape function's derivative of the shape function operation to calculate at Gauss points.
         * 
         */
        void calc_B() 
        {
            for (int i = 0; i < gauss_points_x.size(); ++i)
            {
                B[i](0,0) = -1/length;
                B[i](0,3) = 1/length;
                B[i](1,1) = -6*std::pow(1/length,2) + 12*gauss_points_x[i]*std::pow(1/length,3);
                B[i](1,2) = - 4/length + 6*gauss_points_x[i]*std::pow(1/length, 2);
                B[i](1,4) = 6*std::pow(1/length, 2) - 12*gauss_points_x[i]*std::pow(1/length, 3);
                B[i](1,5) = -2/length + 6 * gauss_points_x[i]* std::pow(1/length,2);
            }
        }

        /**
         * @brief calculates the transformation matrix from the orientation object.
         * 
         * @param sec_offset is the y-axis offset of the section nodes from their actual coordinates.
         * @param origin_x the x-axis of the global coordinate system.
         */
        void calc_T(real sec_offset = 0.0, coords origin_x = {1.0, 0.0, 0.0}) {
            orient.evaluate(nodes, sec_offset, origin_x);   
        }

        /**
         * @brief calculates local constitutive matrix from section information. For this lineaer element it is simply EA and EI along diagonals 1,1 and 2,2
         * 
         */
        void calc_local_constitutive_mat() {
            // given all constitutive mat elements are zeroed we only need to calculate the non-zero diagonal members of this element.
            local_constitutive_mat[0](0,0) = section[0]->get_E()*section[0]->get_A();
            local_constitutive_mat[0](1,1) = section[0]->get_E()*section[0]->get_I();
        }

        /**
         * @brief calcualtes the local strains from the relationship \f$\boldsymbol{\sigma} = \boldsymbol{B}\boldsymbol{d}\f$.
         * @warning requires \ref B to be calculated before this function is called.
         */
        void calc_eps() {
            local_eps[0] = B[0] * local_d;
        }

        /**
         * @brief calculates the local stresses from \f$\boldsymbol{\sigma}=\boldsymbol{D}{\boldsymbol{\varepsilon}}\f$
         * @warning requires \ref local_constitutive_mat and \ref local_eps to be calculated before this function is called.
         */
        void calc_stresses()  {local_stresses[0] = local_constitutive_mat[0]*local_eps[0];}
        
        /**
         * @brief calculates element nodal forces based on nodal displacements and element stiffness.
         * @details calculates the nodal forces from the relationship \f$\boldsymbol{f} = \boldsymbol{k}\boldsymbol{d}\f$
         */
        void calc_local_f() {
            local_f = local_tangent_stiffness*local_d;
        }
                
        /**
         * @brief updates element nodal displacements, strains, stresses, element resistance forces.
         * @warning this is a purely linear element so no re-evaluation of material strength or constitutive matrix.
         * @warning This calculation is NOT done for Gauss points, but only for midpoint of element (which is preferred for output postprocessing).
         */
        void update_state() 
        {
            get_U_from_nodes();
            // need to retrieve local displacement from the global displacement first of all.
            calc_d_from_U();
            // calculating element strain and stress states depends on local displacement d, even though B calculation currently does not.
            calc_B();
            calc_eps();
            calc_stresses();
            // stiffness calculation must come AFTER stress calculation as stiffness may depend on stress.
            calc_stiffnesses();
            // nodal forces depend on stiffness, so must come after stiffness calculations.
            calc_local_f();
            // these local element nodal forces are then transformed to global nodal forces.
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
         * @brief calculates the material stiffness matrix using the standard Euler-Bernouli definition.
         * 
         */
        void calc_mat_stiffness() {
                real A = section[0]->get_A();
                real E = section[0]->get_E();
                real I = section[0]->get_I();
                real EA = E*A;
                real EI = E*I;
                // Row 1
                local_mat_stiffness(0,0) = EA/length;
                local_mat_stiffness(0,3) = -EA/length;
                // Row 2
                local_mat_stiffness(1,1) = 12*EI/std::pow(length,3);
                local_mat_stiffness(1,2) = 6*EI/std::pow(length,2);
                local_mat_stiffness(1,4) = -12*EI/std::pow(length,3);
                local_mat_stiffness(1,5) = 6*EI/std::pow(length,2);
                // Row 3
                local_mat_stiffness(2,1) = 6*EI/std::pow(length,2);
                local_mat_stiffness(2,2) = 4*EI/length;
                local_mat_stiffness(2,4) = -6*EI/std::pow(length,2);
                local_mat_stiffness(2,5) = 2*EI/length;
                // Row 4
                local_mat_stiffness(3,0) = -EA/length;
                local_mat_stiffness(3,3) = EA/length;
                // Row 5
                local_mat_stiffness(4,1) = -12*EI/std::pow(length,3);
                local_mat_stiffness(4,2) = -6*EI/std::pow(length,2);
                local_mat_stiffness(4,4) = 12*EI/std::pow(length,3);
                local_mat_stiffness(4,5) = -6*EI/std::pow(length,2);
                // Row 6
                local_mat_stiffness(5,1) = 6*EI/std::pow(length,2);
                local_mat_stiffness(5,2) = 2*EI/length;
                local_mat_stiffness(5,4) = -6*EI/std::pow(length,2);
                local_mat_stiffness(5,5) = 4*EI/length;
        }

        /**
         * @brief since this is a purely linear element, there is no geometric stiffness to calculate - do nothing.
         * 
         */
        void calc_geom_stiffness() 
        {
            // does nothing. Geometric stiffness is zero.
        }
        /**
         * @brief since this is a purely linear element, there is no external geometric stiffness to calculate - do nothing.
         * 
         */
        virtual void calc_external_geom_stiffness()
        {
            // does nothing. External geometric stiffness is zero.
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
            elem_global_stiffness = orient.get_T().transpose() * local_tangent_stiffness * orient.get_T();
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
         * @brief maps global freedoms to element local freedoms using the transformation matrix.
         * @details uses the relationship \f$\boldsymbol{d} = \boldsymbol{T}\boldsymbol{U}\f$. U comes from nodal displacements.
         * 
         */
        void calc_d_from_U() {local_d = orient.get_T()*global_ele_U;}

        /**
         * @brief Calculates the resistance forces from the relationship \f$ \boldsymbol{R} = \boldsymbol{T}^T\boldsymbol{f}\f$.
         */
        void calc_element_global_resistance_forces() {
            element_global_resistance_forces = orient.get_T().transpose()*local_f;
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
