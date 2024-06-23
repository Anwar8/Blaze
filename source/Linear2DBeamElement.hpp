/**
 * @file Linear2DBeamElement.hpp
 * @brief definitions for the basic Euler-Bernouli linear beam-column.
 */

#ifndef LINEAR_2D_BEAM_ELEMENT_HPP
#define LINEAR_2D_BEAM_ELEMENT_HPP

#include <string>
#include <array>
#include <memory>
#include "main.hpp"
#include "BeamElement2DTemplate.hpp"

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
class Linear2DBeamElement : public BeamElement2DTemplate {
    private:
    protected:
        /**
         * @name basic_information
         * @brief the basic data about the generic beam element
         */
        //@{      
        unsigned id = 0; /**< unique identifier for the element.*/
        std::string const elem_type = "2D_EulerBernouli_beam-column"; /**< string that represents the type of the element.*/
        int ndofs = 3; /**< number of freedoms at each node.*/
        int nnodes = 2; /**< number of nodes.*/
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
         */
        //@{
        vec global_ele_U= make_xd_vec(12); /**< global nodal displacements for all freedoms of the nodes corresponding to the element.*/
        vec local_d = make_xd_vec(6); /**< local nodal-displacements for all freedoms.*/
        vec local_f = make_xd_vec(6); /**< local nodal-forces corresponding to all freedoms.*/
        vec element_resistance_forces = make_xd_vec(12); /**< transformed resistance forces of the element from \ref local_f.*/
        std::vector<spnz> global_R_triplets; /**< triplet vector for global resistance forces \f$\boldsymbol{R}\f$.*/
        vec local_eps = make_xd_vec(2); /**< local strains.*/
        vec local_stresses = make_xd_vec(2); /**< local stresses.*/
        mat local_constitutive_mat = make_xd_mat(2,2); /**< local constitutive matrix \f$\boldsymbol{D}\f$.*/
        mat local_mat_stiffness = make_xd_mat(6,6); /**< local element material stiffness matrix.*/
        mat local_geom_stiffness = make_xd_mat(6,6); /**< local element geometric stiffness matrix.*/
        mat local_tangent_stiffness = make_xd_mat(6,6); /**< local element tangent stiffness matrix.*/
        mat elem_global_stiffness = make_xd_mat(12,12); /**< the global contribution of the element - as in, tangent stiffness after transform via \f$ \boldsymbol{K}_t^e = \boldsymbol{T}^T \boldsymbol{k}_t \boldsymbol{T}\f$*/
        std::vector<spnz> global_stiffness_triplets; /**< the global contributions of the element to the global stiffness - made as sparse matrix contributions that would be gatehred to create the global sparse matrix.*/
        //@}
        

        /**
         * @brief maps local stiffness contributions to their global positions in the stiffness matrix
         * 
         * @details a vector of std array of size 4. the first two indices of the array refer to the
         * transformed local stiffness matrix indices, and the last two refer to the indices where that
         * local stiffness would go in the global stiffness matrix. The size of the std vector is dependent on the
         * number of nodes and which DOFs are active/not fixed. See \ref map_stiffness for details. The size 4 does
         * not change regardless of element type and implementation.
         */
        std::vector<std::array<int,4>> stiffness_map; 
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
            calc_length();
            calc_T();
            calc_local_constitutive_mat();
            calc_stiffnesses();
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
        
        // virtual void calc_N(real x) const = 0;
        /**
         * @brief call the shape function's derivative of the shape function operation to calculate at a specific point.
         * 
         * @param x location along beam-column element at which to calculate the derivative of the shape function.
         */
        void calc_B(real x) {
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
        void calc_nodal_forces() {local_f = local_tangent_stiffness*local_d;}
                
        /**
         * @brief updates element nodal displacements, strains, stresses, element resistance forces.
         * @todo no need to evaluate EA and EI every state update - move to initialisation stage.
         * @todo Change the way B matrix is calculated so it is not done for midpoint of beam.
         * @todo add calculation of geometric and tangent stiffnesses.
         * @todo remove redundant calculation of material stiffness unless it needs to be recalculated again. Remember, it was calculated before finding global U to begin with.
         * @warning calculates \f$\boldsymbol{B}\f$ based on mid-length of the beam not Gauss points.
         */
        virtual void update_state() const = 0;

        
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
        // /**
        //  * @brief runs all stiffness function calculations
        //  * 
        //  */
        // void calc_stiffnesses()
        // {
        //     calc_mat_stiffness();
        //     calc_geom_stiffness();
        //     calc_tangent_stiffness();
        //     calc_elem_global_stiffness();
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
 
        // /**
        //  * @brief prints the internal state of the element.
        //  * 
        //  */
        // void print_element_state(bool print_stresses = true, bool print_strains = false,
        //                          bool print_nodal_disp = false, bool print_nodal_forces = false) 
        // {
        //     if (print_stresses) 
        //     {
        //         std::cout << "element " << id << " stresses are:"  << std::endl << local_stresses << std::endl;
        //     }
        //     if (print_strains) 
        //     {
        //         std::cout << "element " << id << " strains are:"  << std::endl << local_eps << std::endl;
        //     }
        //     if (print_nodal_forces) 
        //     {
        //         std::cout << "element " << id << " nodal forces are:"  << std::endl << local_f << std::endl;
        //     }
        //     if (print_nodal_disp) 
        //     {
        //         std::cout << "element " << id << " nodal displacements are:" << std::endl << local_d << std::endl;
        //     }
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
        void get_U_from_nodes() {
            std::array<real, 6> nodal_disp;
            int i = 0;
            // global_ele_U
            for (auto node: nodes)
            {
                nodal_disp = node->get_nodal_displacements();
                for (auto dof: nodal_disp)
                {
                    global_ele_U(i) = dof;
                    ++i;
                }
            }
        }

        /**
         * @brief Calculates the resistance forces from the relationship \f$ \boldsymbol{R} = \boldsymbol{T}^T\boldsymbol{f}\f$.
         */
        void calc_global_resistance_forces() {
            element_resistance_forces = orient.get_T().transpose()*local_f;
        }

        /**
         * @brief Populates the resistance forces triplets removing any inactive freedoms.
         * @todo I seem to be doing the active_dofs thing too often. The element should also have a set that contains its active dofs!
         */
        virtual void populate_resistance_force_triplets() const = 0;

        /**
         * @brief calculates the global stiffness contribution of the local element and populates global_stiffness_triplets
         * 
         * @details first, the freedoms are mapped to the right size by pre- and post-multiplying by the T matrix
         * After that, \ref stiffness_map is used to map where these contributions would go in the global stiffness
         * matrix. So, this function will populate \ref global_stiffness_triplets with sparse matrix notation
         * 
         */
        virtual void calc_K_global() const = 0;

        /**
         * @brief populates \ref stiffness_map considering active and inactive DOFs for each node of the element
         * 
         * @details see function \ref calc_K_global, and variables \ref stiffness_map, and \ref global_stiffness_triplets. 
         * 
         * @todo REALLY needs to be revisited. attempt to rewrite this function so it does the following:
         *  1. gets all the contribution without worrying about active or not
         *  2. if a contribution is inactive then that contribution is zeroed AND
         *  3. zeroed contributions are not added to \ref global_stiffness_triplets
         * 
         */
        virtual void map_stiffness() const = 0;
        
        /**
         * @brief populate \ref global_dof_map; appears deprecated.
         * 
         */
        virtual void create_dof_map() const = 0;

        /**
         * @brief a function to take care of correctly mapping only active DOFs; appears to have been deprecated.
         * 
         * @param elem_dofs 
         * @param active_dofs 
         * @return std::vector<int> 
         */
        virtual std::vector<int> map_dofs(std::vector<int> elem_dofs, std::set<int> active_dofs) const = 0;
    //@}

    /**
     * @name getter functions
     * @brief functions that retrieve protected variables
     */
    //@{

        int get_ndofs() {return ndofs;}
        mat get_N() {return shape_func.get_N();}
        mat get_B() {return shape_func.get_B();}
        mat get_k() {return shape_func.get_k();}
        mat get_T() {return orient.get_T();}
        vec get_eps() {return local_eps;}
        vec get_d() {return local_d;}

        virtual std::vector<spnz> get_global_resistance_force_triplets() const = 0;
        virtual std::vector<spnz> get_K_global() const = 0;
        virtual int const get_nth_node_id(int n) const = 0;
    //@}

        virtual void set_d(vec new_disp) const = 0;
};


#endif
