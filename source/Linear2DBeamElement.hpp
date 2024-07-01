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
 * 
 */
class Linear2DBeamElement : public BeamElementCommonInterface {
    private:
    protected:
        /**
         * @name basic_information
         * @brief the basic data about the generic beam element
         */
        //@{      
        // unsigned id = 0; /**< unique identifier for the element.*/
        // std::string const elem_type = "2D_EulerBernouli_beam-column"; /**< string that represents the type of the element.*/
        // int ndofs = 3; /**< number of freedoms at each node.*/
        // int nnodes = 2; /**< number of nodes.*/
        // std::vector<real> gauss_points; /**< length-wise coordinates of the Gauss Points. to be set by \ref set_gauss_points*/
        // real length = 0.0; /**< the length for the beam-column element - to be calculated by the orientation object.*/
        //@}

        /**
         * @name beam_basic_objects
         * @brief basic objects needed by the beam-column elements. Section, shape function, transformation, etc.
         */
        //@{
        // std::vector<std::shared_ptr<Node>> nodes; /**< a std::vector that holds the shared ptrs to the nodes.*/
        // BasicSection section; /**< the section for the beam-column element.*/
        // BasicShapeFunction shape_func; /**< the shape function for the beam-column element.*/
        // BasicOrientation orient; /**< the orientation object for the beam-column element.*/
        // NonlinearTransform transformation; /**< the nonlinear transformation used to account for geometric nonlinearity*/
        //@}

        /**
         * @name beam_state_containers
         * @brief the containers for the beam state such as displacement, strain, force, etc.
         * @warning remove the initialisation of the matrices from here and move it to the \ref initialise function!
         */
        // //@{
        // vec global_ele_U= make_xd_vec(12); /**< global nodal displacements for all freedoms of the nodes corresponding to the element.*/
        // vec local_d = make_xd_vec(6); /**< local nodal-displacements for all freedoms.*/
        // vec local_f = make_xd_vec(6); /**< local nodal-forces corresponding to all freedoms.*/
        // vec element_resistance_forces = make_xd_vec(12); /**< transformed resistance forces of the element from \ref local_f.*/
        // std::vector<spnz> global_R_triplets; /**< triplet vector for global resistance forces \f$\boldsymbol{R}\f$.*/
        // vec local_eps = make_xd_vec(2); /**< local strains.*/
        // vec local_stresses = make_xd_vec(2); /**< local stresses.*/
        // mat local_constitutive_mat = make_xd_mat(2,2); /**< local constitutive matrix \f$\boldsymbol{D}\f$.*/
        // mat local_mat_stiffness = make_xd_mat(6,6); /**< local element material stiffness matrix.*/
        // mat local_geom_stiffness = make_xd_mat(6,6); /**< local element geometric stiffness matrix.*/
        // mat local_tangent_stiffness = make_xd_mat(6,6); /**< local element tangent stiffness matrix.*/
        // mat elem_global_stiffness = make_xd_mat(12,12); /**< the global contribution of the element - as in, tangent stiffness after transform via \f$ \boldsymbol{K}_t^e = \boldsymbol{T}^T \boldsymbol{k}_t \boldsymbol{T}\f$*/
        // std::vector<spnz> global_stiffness_triplets; /**< the global contributions of the element to the global stiffness - made as sparse matrix contributions that would be gatehred to create the global sparse matrix.*/
        // // @}
        

        /**
         * @brief maps local stiffness contributions to their global positions in the stiffness matrix
         * 
         * @details a vector of std array of size 4. the first two indices of the array refer to the
         * transformed local stiffness matrix indices, and the last two refer to the indices where that
         * local stiffness would go in the global stiffness matrix. The size of the std vector is dependent on the
         * number of nodes and which DOFs are active/not fixed. See \ref map_stiffness for details. The size 4 does
         * not change regardless of element type and implementation.
         */
        // std::vector<std::array<int,4>> stiffness_map; 
        //@}

        

    public:
    /**
     * @name element initialisers
     * @brief functions that deal with constructing and initialising the beam-column element
     */
    //@{
        /**
         * @brief Construct a new 2D Beam Element object.
         * 
         * @tparam Container any type of std container that has a std::size and built-in iterators
         * @param given_id unique identifier for the element; will be passed to the nodes
         * @param in_nodes a container of shared pointers to node objects
         */
        template<typename Container>
        Linear2DBeamElement(int given_id, Container& in_nodes) {
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
            elem_type = "2D_EulerBernouli_beam-column"; /**< string that represents the type of the element.*/
            ndofs = 3; /**< number of freedoms at each node. 3 for this element type.*/
            nnodes = 2; /**< number of nodes. 2 nodes for this element type.*/
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
            for (auto node : in_nodes) {
                
                node->add_connected_element(id);
            }
            calc_T();
            calc_length();
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
            element_resistance_forces = make_xd_vec(12); /**< transformed resistance forces of the element from \ref local_f. 6 per node since we have 3D global nodes, and two nodes giving us 12 components.*/
            local_eps = make_xd_vec(2); /**< local strains. \f$ \boldsymbol{\varepsilon} = \begin{bmatrix} \varepsilon_{xx} & \kappa\end{bmatrix}^T\f$*/
            local_stresses = make_xd_vec(2); /**< local stresses.\f$ \boldsymbol{\sigma} = \begin{bmatrix} N & M \end{bmatrix}^T\f$*/
            local_constitutive_mat = make_xd_mat(2,2); /**< local constitutive matrix \f$\boldsymbol{D} = \begin{bmatrix} EA & 0 \\ 0 & EI\end{bmatrix}\f$.*/
            local_mat_stiffness = make_xd_mat(6,6); /**< local element material stiffness matrix.*/
            local_geom_stiffness = make_xd_mat(6,6); /**< local element geometric stiffness matrix.*/
            local_tangent_stiffness = make_xd_mat(6,6); /**< local element tangent stiffness matrix.*/
            elem_global_stiffness = make_xd_mat(12,12); /**< the global contribution of the element - as in, tangent stiffness after transform via \f$ \boldsymbol{K}_t^e = \boldsymbol{T}^T \boldsymbol{k}_t \boldsymbol{T}\f$*/
        }

        /**
         * @brief Set the gauss points std::vector.
         * @details Sets \ref gauss_points to an appropriate size and set of values. To be called by constructor.
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
         * @brief retrieves the length of the element using the orientation object since this is a linear element.
         * 
         */
        void calc_length() {
            length = orient.get_length();
        }
        /**
         * @brief calculates the shape function of the beam-column element for location x.
         * 
         * @param x location at which shape function is evaluated.
         */
        void calc_N(real x)
        {
            shape_func.calc_N(x, length);
        }
        /**
         * @brief call the shape function's derivative of the shape function operation to calculate at a specific point.
         * 
         * @param x location along beam-column element at which to calculate the derivative of the shape function.
         */
        void calc_B(real x) 
        {
            shape_func.calc_B(x, length);
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
         * @brief calculates local constitutive matrix from section information.
         * 
         */
        void calc_local_constitutive_mat() {
            // given all constitutive mat elements are zeroed we only need to calculate the non-zero diagonal members of this element.
            local_constitutive_mat(0,0) = section.get_E()*section.get_A();
            local_constitutive_mat(1,1) = section.get_E()*section.get_I();
        }

        /**
         * @brief calcualtes the local strains from the relationship \f$\boldsymbol{\sigma} = \boldsymbol{B}\boldsymbol{d}\f$.
         * 
         */
        void calc_eps() {
            local_eps = shape_func.get_B() * local_d;
        }

        /**
         * @brief calculates the local stresses from \f$\boldsymbol{\sigma}=\boldsymbol{D}{\boldsymbol{\varepsilon}}\f$
         * @warning depends on `Eigen3` overlay for the \* operation for matrix objects. 
         */
        void calc_stresses()  {local_stresses = local_constitutive_mat*local_eps;}
        
        /**
         * @brief calculates element nodal forces based on nodal displacements and element stiffness.
         * @details calculates the nodal forces from the relationship \f$\boldsymbol{f} = \boldsymbol{k}\boldsymbol{d}\f$
         */
        void calc_local_f() {local_f = local_tangent_stiffness*local_d;}
                
        /**
         * @brief updates element nodal displacements, strains, stresses, element resistance forces.
         * @warning this is a purely linear element so no re-evaluation of material strength or constitutive matrix.
         * @warning This calculation is NOT done for Gauss points, but only for midpoint of element (which is preferred for output postprocessing).
         */
        void update_state() 
        {
            // need to retrieve local displacement from the global displacement first of all.
            calc_d_from_U();
            // calculating element strain and stress states depends on local displacement d, even though B calculation currently does not.
            calc_B(length*0.5);
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
         * @todo double check whether this is the Euler-Bernouli definition or a Timoshenko definition.
         * 
         */
        void calc_mat_stiffness() {
                real A = section.get_A();
                real E = section.get_E();
                real I = section.get_I();
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
        /**
         * @brief runs all stiffness function calculations
         * 
         */
        // void calc_stiffnesses()
        // {
        //     this->BeamElementBaseClass::calc_stiffnesses();
        // }
    //@}

    /**
     * @name logging functions
     * @brief functions used for logging output to stream - mostly for debugging.
     */
    //@{
        /**
         * @brief prints the most important information of the element to the output stream.
         */
        void print_info() {
            std::cout << "elem " << id << " of type " <<elem_type << " with " << ndofs << " dofs, and " << nnodes << " nodes:" << std::endl;
            for (auto node_i: nodes) {
                node_i->print_info();
            }
            std::cout << "it is also of length " << length << std::endl;
        }
 
        /**
         * @brief prints the internal state of the element.
         * 
         */
        // void print_element_state(bool print_stresses = true, bool print_strains = false,
        //                          bool print_nodal_disp = false, bool print_nodal_forces = false) 
        // {
        //   this->BeamElementBaseClass::print_element_state(print_stresses, print_strains, print_nodal_disp, print_nodal_forces);
        // }
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
         * @brief Get the \ref global_ele_U from each node object connected to the element.
         * 
         */
        // void get_U_from_nodes() 
        // {
        //     this->BeamElementBaseClass::get_U_from_nodes();
        // }

        /**
         * @brief Calculates the resistance forces from the relationship \f$ \boldsymbol{R} = \boldsymbol{T}^T\boldsymbol{f}\f$.
         */
        void calc_element_global_resistance_forces() {
            element_resistance_forces = orient.get_T().transpose()*local_f;
        }

        /**
         * @brief implemented in base class.
         */
        // void populate_resistance_force_triplets() {
        //     this->BeamElementBaseClass::populate_resistance_force_triplets();
        // }

        /**
         * @brief calculates the global stiffness contribution of the local element and populates global_stiffness_triplets
         * 
         * @details first, the freedoms are mapped to the right size by pre- and post-multiplying by the T matrix
         * After that, \ref stiffness_map is used to map where these contributions would go in the global stiffness
         * matrix. So, this function will populate \ref global_stiffness_triplets with sparse matrix notation
         * 
         */
        // void calc_global_stiffness_triplets() {
        //     this->BeamElementBaseClass::calc_global_stiffness_triplets();
        // }

        /**
         * @brief populates \ref stiffness_map considering active and inactive DOFs for each node of the element
         * 
         * @details see function \ref calc_global_stiffness_triplets, and variables \ref stiffness_map, and \ref global_stiffness_triplets. 
         * 
         * @todo REALLY needs to be revisited. attempt to rewrite this function so it does the following:
         *  1. gets all the contribution without worrying about active or not
         *  2. if a contribution is inactive then that contribution is zeroed AND
         *  3. zeroed contributions are not added to \ref global_stiffness_triplets
         * 
         */
        // void map_stiffness() 
        // {
        //     this->BeamElementBaseClass::map_stiffness();
        // }

        /**
         * @brief a function to take care of correctly mapping only active DOFs; appears to have been deprecated.
         * 
         * @param elem_dofs 
         * @param active_dofs 
         * @return std::vector<int> 
         */
        // std::vector<int> map_dofs(std::vector<int> elem_dofs, std::set<int> active_dofs)
        // {
        //     return this->BeamElementBaseClass::map_dofs(elem_dofs, active_dofs);
        // }
    //@}

    /**
     * @name getter functions
     * @brief functions that retrieve protected variables
     */
    //@{

    //     int get_ndofs() {return this->BeamElementBaseClass::get_ndofs();}
    //     mat get_N() {return this->BeamElementBaseClass::get_N();}
    //     mat get_B() {return this->BeamElementBaseClass::get_B();}
    //     mat get_k() {return this->BeamElementBaseClass::get_k();}
    //     mat get_T() {return this->BeamElementBaseClass::get_T();}
    //     vec get_eps() {return this->BeamElementBaseClass::get_eps();}
    //     vec get_d() {return this->BeamElementBaseClass::get_d();}

    //     std::vector<spnz> get_global_resistance_force_triplets() {return this->BeamElementBaseClass::get_global_resistance_force_triplets();}
    //     virtual std::vector<spnz> get_K_global() {return this->BeamElementBaseClass::get_K_global();}
    //     int const get_nth_node_id(int n) {return this->BeamElementBaseClass::get_nth_node_id(n);}
    // //@}

    //     void set_d(vec new_disp) {this->BeamElementBaseClass::set_d(new_disp);}
};


#endif
