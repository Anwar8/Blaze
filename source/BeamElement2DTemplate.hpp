/**
 * @file beam_element.hpp
 * @brief definitions for basic beam-column elements.
 * @todo get rid of the shape function object, or modify it so that it is more powerful.
 */

#ifndef BEAM_ELEMENT_TEMPLATE_HPP
#define BEAM_ELEMENT_TEMPLATE_HPP

#include <string>
#include <array>
#include <memory>
#include "main.hpp"
#include "maths_defaults.hpp"
#include "basic_section.hpp"
#include "basic_shape_function.hpp"
#include "basic_orientation.hpp"
#include "node.hpp"

/**
 * @brief a 2D beam element base class with unknown total freedoms
 * 
 * @details
 * This base class does not know how many freedoms the element will have, and will thus not know the size of any of the matrices.
 * 
 * @todo 1. add functionality to apply distributed mechanical loading
 *       2. add functionality to apply thermal gradient and uniform temperature
 *       3. add functionality to apply non-uniform thermal loading interpolation (LOW PRIORITY)
 * 
 */
class BeamElement2DTemplate {
    private:
    protected:
        /**
         * @name basic_information
         * @brief the basic data about the generic beam element
         */
        //@{      
        unsigned id = 0; /**< unique identifier for the element.*/
        std::string const elem_type = "beam-column"; /**< string that represents the type of the element.*/
        int ndofs = -1; /**< number of freedoms at each node. Unknown in this virtual class. Set to -1 to force implementation in subclass.*/
        int nnodes = -1; /**< number of nodes. Set to -1 to force implementation in subclass.*/
        std::vector<real> gauss_points; /**< length-wise coordinates of the Gauss Points. to be set by \ref set_gauss_points*/
        real length = 0.0; /**< the length for the beam-column element - to be calculated by the orientation object.*/
        //@}

        /**
         * @name beam_basic_objects
         * @brief basic objects needed by the beam-column elements. Section, shape function, transformation, etc.
         */
        //@{
        std::vector<std::shared_ptr<Node>> nodes; /**< a std::vector that holds the shared ptrs to the nodes.*/
        BasicSection section; /**< the section for the beam-column element.*/
        BasicShapeFunction shape_func; /**< the shape function for the beam-column element.*/
        BasicOrientation orient; /**< the orientation object for the beam-column element.*/
        //@}

        /**
         * @name beam_state_containers
         * @brief the containers for the beam state such as displacement, strain, force, etc.
         */
        //@{
        vec global_ele_U; /**< global nodal displacements for all freedoms of the nodes corresponding to the element.*/
        vec local_d; /**< local nodal-displacements for all freedoms.*/
        vec local_f; /**< local nodal-forces corresponding to all freedoms.*/
        vec element_resistance_forces; /**< transformed resistance forces of the element from \ref local_f.*/
        std::vector<spnz> global_R_triplets; /**< triplet vector for global resistance forces \f$\boldsymbol{R}\f$.*/
        vec local_eps; /**< local strains.*/
        vec local_stresses; /**< local stresses.*/
        mat local_constitutive_mat; /**< local constitutive matrix \f$\boldsymbol{D}\f$.*/
        mat local_mat_stiffness; /**< local element material stiffness matrix.*/
        mat local_geom_stiffness; /**< local element geometric stiffness matrix.*/
        mat local_tangent_stiffness; /**< local element tangent stiffness matrix.*/
        mat elem_global_stiffness; /**< the global contribution of the element - as in, tangent stiffness after transform via \f$ \boldsymbol{K}_t^e = \boldsymbol{T}^T \boldsymbol{k}_t \boldsymbol{T}\f$*/
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
         * @brief Construct a new 2D Beam Element object.
         * 
         * @tparam Container any type of std container that has a std::size and built-in iterators
         * @param given_id unique identifier for the element; will be passed to the nodes
         * @param in_nodes a container of shared pointers to node objects
         */
        template<typename Container>
        BeamElement2DTemplate(int given_id, Container& in_nodes) {
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
            std::cout << "This is a pseudo virtual function. Must define an overriding initialiser for your class." << std::endl;
            std::exit(1);  
        }

        /**
         * @brief Set the gauss points std::vector.
         * @details Sets \ref gauss_points to an appropriate size and set of values. To be called by constructor.
         * 
         */
        virtual void set_gauss_points() const = 0;
        
        /**
         * @brief prints the most important information of the element to the output stream.
         */
        virtual void print_info() const = 0;
        /**
         * @brief 
         * 
         */
        virtual void calc_length() const = 0;

        // virtual void calc_N(real x) const = 0;
        /**
         * @brief call the shape function's derivative of the shape function operation to calculate at a specific point.
         * 
         * @param x 
         */
        void calc_B(real x);
        /**
         * @brief calculates the material stiffness matrix using the shape-function function \ref ShapeFunction::calc_elem_mat_stiffness.
         * 
         */
        void calc_mat_stiffness() {shape_func.calc_elem_mat_stiffness(length, section, local_mat_stiffness);}
        /**
         * @brief calculates the geometric stiffness matrix using the shape-function function \ref ShapeFunction::calc_elem_geom_stiffness.
         * 
         */
        void calc_geom_stiffness() {shape_func.calc_elem_geom_stiffness(length, local_stresses(0), local_geom_stiffness);}
        /**
         * @brief sums the \ref local_mat_stiffness and \ref local_geom_stiffness to create the tangent stiffness.
         * 
         */
        void calc_tangent_stiffness() {
            std::cout << ">>>>>>>>>>>calc_tangent_stifness<<<<<<<<<<" << std::endl;
            std::cout << "calc_tangent_stiffness::elem " << id << " mat_stiffness is " << std::endl << local_mat_stiffness << std::endl;
            std::cout << "calc_tangent_stiffness::elem " << id << " geom_stiffness is " << std::endl << local_geom_stiffness << std::endl;
            std::cout << "calc_tangent_stiffness::elem " << id << " mat_stiffness + geom_stiffness is " << std::endl << local_mat_stiffness + local_geom_stiffness << std::endl;
            this->local_tangent_stiffness = this->local_mat_stiffness + this->local_geom_stiffness;
            std::cout << "calc_tangent_stiffness::elem " << id << " tangent_stiffness is " << std::endl << local_tangent_stiffness << std::endl;
        }

        /**
         * @brief calculates the global contributions of the tangent stiffness to the global stiffness matrix as a 6x6 matrix.
         * 
         */
        void calc_elem_global_stiffness() {
            elem_global_stiffness = orient.get_T().transpose() * local_tangent_stiffness * orient.get_T();
        }
        /**
         * @brief calculates the transformation matrix from the orientation object.
         * 
         * @param sec_offset is the y-axis offset of the section nodes from their actual coordinates.
         * @param origin_x the x-axis of the global coordinate system.
         */
        void calc_T(real sec_offset = 0.0, coords origin_x = {1.0, 0.0, 0.0});

        /**
         * @brief calculates local constitutive matrix from section information.
         * 
         */
        void calc_local_constitutive_mat();

        /**
         * @brief calcualtes the local strains from the relationship \f$\boldsymbol{\sigma} = \boldsymbol{B}\boldsymbol{d}\f$.
         * 
         */
        void calc_eps(){
            local_eps = shape_func.get_B() * local_d;
            if (VERBOSE)
            {
                std::cout << "Element " << id << " B is " << std::endl << shape_func.get_B() << std::endl;
                std::cout << "Element " << id << " d is " << std::endl << local_d << std::endl;
            }
            }
        /**
         * @brief calculates the local stresses from \f$\boldsymbol{\sigma}=\boldsymbol{D}{\boldsymbol{\varepsilon}}\f$
         * @warning depends on `Eigen3` overlay for the \* operation for matrix objects. 
         */
        void calc_stresses() {local_stresses = local_constitutive_mat*local_eps;}
        
        /**
         * @brief calculates element nodal forces based on nodal displacements and element stiffness.
         * @details calculates the nodal forces from the relationship \f$\boldsymbol{f} = \boldsymbol{k}\boldsymbol{d}\f$
         */
        void calc_nodal_forces() {local_f = local_tangent_stiffness*local_d;}

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
        void get_U_from_nodes();
        
        /**
         * @brief updates element nodal displacements, strains, stresses, element resistance forces.
         * @todo no need to evaluate EA and EI every state update - move to initialisation stage.
         * @todo Change the way B matrix is calculated so it is not done for midpoint of beam.
         * @todo add calculation of geometric and tangent stiffnesses.
         * @todo remove redundant calculation of material stiffness unless it needs to be recalculated again. Remember, it was calculated before finding global U to begin with.
         * @warning calculates \f$\boldsymbol{B}\f$ based on mid-length of the beam not Gauss points.
         */
        void update_state()
        {

         // really should not be doing this for midpoint location.
         calc_B(length*0.5);

         get_U_from_nodes();
         calc_d_from_U();
         calc_eps();
         calc_stresses();
         calc_mat_stiffness();
         calc_geom_stiffness();
         calc_tangent_stiffness();
         calc_elem_global_stiffness();
         calc_nodal_forces();
         populate_resistance_force_triplets();
        }

        void calc_stiffnesses()
        {
            calc_mat_stiffness();
            calc_geom_stiffness();
            calc_tangent_stiffness();
            calc_elem_global_stiffness();
        }

        /**
         * @brief prints the internal state of the element.
         * 
         */
        void print_element_state(bool print_stresses = true, bool print_strains = false,
                                 bool print_nodal_disp = false, bool print_nodal_forces = false) 
        {
            if (print_stresses) 
            {
                std::cout << "element " << id << " stresses are:"  << std::endl << local_stresses << std::endl;
            }
            if (print_strains) 
            {
                std::cout << "element " << id << " strains are:"  << std::endl << local_eps << std::endl;
            }
            if (print_nodal_forces) 
            {
                std::cout << "element " << id << " nodal forces are:"  << std::endl << local_f << std::endl;
            }
            if (print_nodal_disp) 
            {
                std::cout << "element " << id << " nodal displacements are:" << std::endl << local_d << std::endl;
            }
        }
        /**
         * @brief calculates the global stiffness contribution of the local element and populates global_stiffness_triplets
         * 
         * @details first, the freedoms are mapped to the right size by pre- and post-multiplying by the T matrix
         * After that, \ref stiffness_map is used to map where these contributions would go in the global stiffness
         * matrix. So, this function will populate \ref global_stiffness_triplets with sparse matrix notation
         * 
         */
        void calc_K_global();

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
        void map_stiffness();
        
        /**
         * @brief populate \ref global_dof_map; appears deprecated.
         * 
         */
        void create_dof_map();

        /**
         * @brief a function to take care of correctly mapping only active DOFs; appears to have been deprecated.
         * 
         * @param elem_dofs 
         * @param active_dofs 
         * @return std::vector<int> 
         */
        std::vector<int> map_dofs(std::vector<int> elem_dofs, std::set<int> active_dofs);
        

        int get_ndofs() {return ndofs;}
        mat get_N() {return shape_func.get_N();}
        mat get_B() {return shape_func.get_B();}
        mat get_k() {return shape_func.get_k();}
        mat get_T() {return orient.get_T();}
        vec get_eps() {return local_eps;}
        vec get_d() {return local_d;}

        /**
         * @brief Calculates the resistance forces from the relationship \f$ \boldsymbol{R} = \boldsymbol{T}^T\boldsymbol{f}\f$.
         */
        void calc_global_resistance_forces()
        {
            element_resistance_forces = orient.get_T().transpose()*local_f;
        }
        /**
         * @brief Populates the resistance forces triplets removing any inactive freedoms.
         * @todo I seem to be doing the active_dofs thing too often. The element should also have a set that contains its active dofs!
         */
        void populate_resistance_force_triplets() {
            global_R_triplets.clear();
            // the 12x1 full resistance vector from local nodal forces vector f
            if (VERBOSE)
            {
            std::cout << "element " << id << " has untransformed local_f " <<std::endl << local_f << std::endl;
            }
            vec element_resistance_forces = orient.get_T().transpose()*local_f;
            if (VERBOSE)
            {
            std::cout << "element " << id << " has element_resistance_forces " <<std::endl << element_resistance_forces << std::endl;
            }
            std::set<int> node_active_dofs;
            int nz_i = 0;
            real force_value;
            int total_nodal_ndofs_completed = 0; // each node we finish with, we add 6 to this. 
            // This means we have to move to the next set of values corresponding to the next 
            // node in the full resistance vector.         
            for (auto node: nodes)
            {
                int nodal_dof_index = 0;
                node_active_dofs = node->get_active_dofs();
                nz_i = node->get_nz_i();
                for (auto active_dof: node_active_dofs)
                {
                
                    force_value = element_resistance_forces(active_dof + total_nodal_ndofs_completed);
                    // since inactive nodes do not appear in R, we have to make sure to be careful about where we add our nodal forces.
                    // here, nz_i + nodal_dof_index simply starts at where the node freedoms start in the global index, and then
                    // iterates one by one. See how we ++ nodal_dof_index for each freedom we add, and how we restrat from zero when
                    // we start work with the next node?
                    global_R_triplets.push_back(spnz(nz_i + nodal_dof_index, 0, force_value));
                    if (VERBOSE)
                    {
                    std::cout << "element " << id << " node " << node->get_id() << " pushed force_value " << force_value << " to " << nz_i + nodal_dof_index <<  std::endl; 
                    }
                    nodal_dof_index++;
                }
                //**< has to be 6 because each node has 6 dofs and our \ref element_resistance_forces also has 6 rows for each node!*
                total_nodal_ndofs_completed += 6;
            }
        }
        std::vector<spnz> get_global_resistance_force_triplets() {return global_R_triplets;}

        std::vector<int> get_global_dof_map() {return global_dof_map;}
        std::vector<spnz> get_K_global() {return global_stiffness_triplets;}

        int const get_nth_node_id(int n) const;
        /**
         * @brief moves the location of the nodes up, does NOT cause a transformation. 
         * 
         * @attention This is a setter function, nothing more!
         * 
         * @param up distance to move nodes by, which is done along the y-axis.
         */
        void move_nodes_up(real up) {
            for (auto node: nodes) {
                node->set_z(up);
            }
        }
        void set_d(vec new_disp) {local_d = new_disp;}
};


#endif
