/**
 * @file BeamElementBaseClass.hpp
 * @brief definitions for base class for all basic beam-column elements.
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
#include "NonlinearTransform.hpp"
#include "node.hpp"

/**
 * @brief a 2D beam element base class with unknown total freedoms
 * 
 * @details
 * This base class does not know how many freedoms the element will have, and will thus not know the size of any of the matrices.
 * This is the most fundamental beam-column element class and is inherited by all other beam-column elements implicitly by being the parent of \ref BeamElementCommonInterface.
 * 
 * @todo 1. add functionality to apply distributed mechanical loading
 * @todo 2. add functionality to apply thermal gradient and uniform temperature
 * @todo 3. add functionality to apply non-uniform thermal loading interpolation (LOW PRIORITY)
 * 
 */
class BeamElementBaseClass {
    private:
    protected:
        /**
         * @name basic_information
         * @brief the basic data about the generic beam element
         */
        //@{      
        unsigned id = 0; /**< unique identifier for the element.*/
        std::string elem_type = "pure-virtual-beam-column"; /**< string that represents the type of the element.*/
        int ndofs = -1; /**< number of freedoms at each node. Unknown in this virtual class. Set to -1 to force implementation in subclass.*/
        int nnodes = -1; /**< number of nodes. Set to -1 to force implementation in subclass.*/
        std::vector<real> gauss_points_x; /**< length-wise coordinates of the Gauss Points. to be set by \ref initialise_gauss_points*/
        std::vector<real> gauss_points_w; /**< Weight of the Gauss Points. to be set by \ref initialise_gauss_points*/
        real length = 0.0; /**< the length for the beam-column element - to be calculated by the orientation object.*/
        //@}

        /**
         * @name beam_basic_objects
         * @brief basic objects needed by the beam-column elements. Section, shape function, transformation, etc.
         */
        //@{
        std::vector<std::shared_ptr<Node>> nodes; /**< a std::vector that holds the shared ptrs to the nodes.*/
        BasicSection section; /**< the section for the beam-column element.*/
        BasicOrientation orient; /**< the orientation object for the beam-column element.*/
        NonlinearTransform transformation; /**< the nonlinear transformation used to account for geometric nonlinearity*/
        //@}

        /**
         * @name beam_state_containers
         * @brief the containers for the beam state such as displacement, strain, force, etc.
         * @todo should the shape function's \f$ \boldysymbol{N}\f$ and \f$ \boldysymbol{B}\f$ be placed here? if so, are they calculated at the Gauss points? 
         */
        //@{
        vec global_ele_U; /**< global nodal displacements for all freedoms of the nodes corresponding to the element.*/
        vec local_d; /**< local nodal-displacements for all freedoms.*/
        vec local_f; /**< local nodal-forces corresponding to all freedoms.*/
        vec element_global_resistance_forces; /**< transformed resistance forces of the element from \ref local_f.*/
        std::vector<spnz> global_R_triplets; /**< triplet vector for global resistance forces \f$\boldsymbol{R}\f$.*/
        std::vector<vec> local_eps; /**< local strains.*/
        std::vector<vec> local_stresses; /**< local stresses.*/
        std::vector<mat> N; /**< the shape function of the element.*/
        std::vector<mat> B; /**< the derivative of the shape function of the element.*/
        std::vector<mat> local_constitutive_mat; /**< local constitutive matrix \f$\boldsymbol{D}\f$.*/
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
         * @param sect the \ref BasicSection object that contains the material properties of the element.
         */
        template<typename Container>
        BeamElementBaseClass(int given_id, Container& in_nodes, BasicSection sect) {
            initialise(given_id, in_nodes);
        }
        /**
         * @brief Construct a new Beam Element Base Class object that is purely virtual. Needed for inheritance
         * 
         */
        BeamElementBaseClass() {};

        /**
         * @brief initialises the beam column-element with an id, nodes, and initialise any parameters needed for a new element.
         * 
         * @tparam Container any type of std container that has a std::size and built-in iterators
         * @param given_id unique identifier for the element; will be passed to the nodes
         * @param in_nodes a container of shared pointers to node objects
         */
        template<typename Container>
        void initialise(int given_id, Container& in_nodes, BasicSection section) {
            std::cout << "This is a pseudo virtual function. Must define an overriding initialiser for your class." << std::endl;
            std::exit(1);  
        }

        /**
         * @brief initialises the beam state containers (matricces and vectors) to zeros of the right size for the element.
         * 
         */
        virtual void initialise_state_containers() = 0;

        /**
         * @brief Set the gauss points std::vector.
         * @details Sets \ref gauss_points_x and \ref gauss_points_w to an appropriate size and initial set of values. To be called by constructor.
         * 
         */
        virtual void initialise_gauss_points() = 0;
        
        /**
         * @brief Updates Gauss points after the length of the element is known.
         * @details Multiplies \ref gauss_points_x and \ref gauss_points_w by the length of the element which is now known.
         * 
         */
        virtual void update_gauss_points() = 0;
    //@}
    /**
     * @name element operator overloads
     * @brief allows for overloading some operators for easier operations and for some common STL algorithms such as for sorting.
     * 
     */
    //@{
    /**
     * @brief overloads the less than operator to compare elements by their ID, allowing easy sorting of element STL containers via \ref std::sort.
     * @todo Understand how I could implement a pure virtual operator override knowing that the function signatures would not match as the derived 
     * class cannot access the base class protected members, and so must use a operator<(const DerivedClass& other_elem) signature vs the base class' operator<(const BaseClass& other_elem) signature.
     */
    // virtual bool operator<(const BeamElementBaseClass& other_elem) const = 0;
    //@}

    /**
     * @name element property calculation functions
     * @brief functions that are used to calculate element properties such as stress, shape, length, and strain
     */
    //@{
    
        /**
         * @brief calculates the length of the element using the orientation or transformation object
         * 
         */
        virtual void calc_length() = 0;
        /**
         * @brief calculates the shape function of the beam-column element for location of Gauss points.
         * 
         */
        virtual void calc_N() = 0;
        /**
         * @brief call the shape function's derivative of the shape function operation to calculate at location of Gauss points.
         * 
         */
        virtual void calc_B() = 0;

        /**
         * @brief calculates the transformation matrix from the orientation object.
         * 
         * @param sec_offset is the y-axis offset of the section nodes from their actual coordinates.
         * @param origin_x the x-axis of the global coordinate system.
         */
        virtual void calc_T(real sec_offset = 0.0, coords origin_x = {1.0, 0.0, 0.0}) = 0;

        /**
         * @brief calculates local constitutive matrix from section information.
         * 
         */
        virtual void calc_local_constitutive_mat() = 0;

        /**
         * @brief calcualtes the local strains from the relationship \f$\boldsymbol{\sigma} = \boldsymbol{B}\boldsymbol{d}\f$.
         * 
         */
        virtual void calc_eps() = 0;
        /**
         * @brief calculates the local stresses from \f$\boldsymbol{\sigma}=\boldsymbol{D}{\boldsymbol{\varepsilon}}\f$
         * @warning depends on `Eigen3` overlay for the \* operation for matrix objects. 
         */
        virtual void calc_stresses() = 0;
        
        /**
         * @brief calculates element nodal forces based on nodal displacements and element stiffness.
         * @details calculates the nodal forces from the relationship \f$\boldsymbol{f} = \boldsymbol{k}\boldsymbol{d}\f$
         */
        virtual void calc_local_f() = 0;
                
        /**
         * @brief updates element nodal displacements, strains, stresses, element resistance forces.
         * @todo no need to evaluate EA and EI every state update - move to initialisation stage.
         * @todo Change the way B matrix is calculated so it is not done for midpoint of beam.
         * @todo add calculation of geometric and tangent stiffnesses.
         * @todo remove redundant calculation of material stiffness unless it needs to be recalculated again. Remember, it was calculated before finding global U to begin with.
         * @warning calculates \f$\boldsymbol{B}\f$ based on mid-length of the beam not Gauss points.
         */
        virtual void update_state() = 0;


    //@}

    /**
     * @name stiffness matrix functions
     * @brief functions that deal with generating and evaluating the different stiffness matrices.
     */
    //@{
        /**
         * @brief calculates the material stiffness matrix using the shape-function function \ref ShapeFunction::calc_elem_mat_stiffness.
         * 
         */
        virtual void calc_mat_stiffness() = 0;
        /**
         * @brief calculates the geometric stiffness matrix using the shape-function function \ref ShapeFunction::calc_elem_geom_stiffness.
         * 
         */
        virtual void calc_geom_stiffness() = 0;
        /**
         * @brief calculates the "external" nonlinear contribution of the element to the global tangent stiffness.
         * 
         */
        virtual void calc_external_geom_stiffness() = 0;
        /**
         * @brief sums the \ref local_mat_stiffness and \ref local_geom_stiffness to create the tangent stiffness.
         * 
         */
        virtual void calc_tangent_stiffness() = 0;

        /**
         * @brief calculates the global contributions of the tangent stiffness to the global stiffness matrix as a 6x6 matrix.
         * 
         */
        virtual void calc_elem_global_stiffness() = 0;
    
        /**
         * @brief runs all stiffness function calculations
         * 
         */
        virtual void calc_stiffnesses() = 0;
    //@}

    /**
     * @name logging functions
     * @brief functions used for logging output to stream - mostly for debugging.
     */
    //@{
        /**
         * @brief prints the most important information of the element to the output stream.
         */
        virtual void print_info() = 0;
 
        /**
         * @brief prints the internal state of the element.
         * 
         */
        virtual void print_element_state(bool print_stresses = true, bool print_strains = false,
                                 bool print_nodal_disp = false, bool print_nodal_forces = false) = 0;
    //@}

    /**
     * @name local-global mapping function
     * @brief functions that deal mapping degrees of freedoms between local element level and global matrices.
     */
    //@{
        /**
         * @brief Get the \ref global_ele_U from each node object connected to the element.
         * 
         */
        virtual void get_U_from_nodes() = 0;

        /**
         * @brief maps global freedoms to element local freedoms using the transformation matrix.
         * @details uses the relationship \f$\boldsymbol{d} = \boldsymbol{T}\boldsymbol{U}\f$. U comes from nodal displacements.
         * 
         */
        virtual void calc_d_from_U() = 0;

        /**
         * @brief Calculates the resistance forces from the relationship \f$ \boldsymbol{R} = \boldsymbol{T}^T\boldsymbol{f}\f$.
         */
        virtual void calc_element_global_resistance_forces() = 0;

        /**
         * @brief Populates the resistance forces triplets removing any inactive freedoms.
         */
        virtual void populate_resistance_force_triplets() = 0;

        /**
         * @brief calculates the global stiffness contribution of the local element and populates global_stiffness_triplets
         * 
         */
        virtual void calc_global_stiffness_triplets() = 0;


        /**
         * @brief populates \ref stiffness_map considering active and inactive DOFs for each node of the element
         * 
         */
        virtual void map_stiffness() = 0;
        
    //@}
    /**
     * @name setter functions
     * @brief functions that set protected variables
     */
    //@{
        /**
         * @brief Set \ref global_ele_U to some value for testing.
         * 
         * @param global_U_vec a vector that the object's \ref global_ele_U will be replaced by.
         */
        virtual void set_global_U(vec global_U_vec) = 0;
        
        /**
         * @brief Set \ref local_d to some displacement vector.
         * 
         * @param new_disp the new displacement the \ref local_d would be replaced by.
         */
        virtual void set_d(vec new_disp) = 0;
    //@}

    /**
     * @name getter functions
     * @brief functions that retrieve protected variables. Also needed for testing.
     */
    //@{

        virtual int get_ndofs() const = 0;
        virtual int get_nnodes() const = 0;
        virtual std::string get_elem_type() const = 0;
        virtual unsigned get_id() const = 0;

        virtual vec get_global_ele_U() const = 0;
        virtual vec get_local_d() const = 0;
        virtual vec get_local_f() const = 0;
        virtual vec get_element_resistance_forces() const = 0;
        virtual std::vector<spnz> get_global_resistance_force_triplets() = 0;
        virtual vec get_eps() const = 0;
        virtual vec get_local_stresses() const = 0;
        virtual mat get_local_constitutive_mat() const = 0;
        virtual mat get_local_mat_stiffness() const = 0;
        virtual mat get_local_geom_stiffness() const = 0;
        virtual mat get_local_tangent_stiffness() const = 0;
        virtual mat get_elem_global_stiffness() const = 0;
        virtual std::vector<spnz> get_global_stiffness_triplets() = 0;

        virtual mat get_N() const = 0;
        virtual mat get_B() const = 0;
        virtual mat get_T() = 0;
        virtual real get_L() = 0;
        
        virtual int const get_nth_node_id(int n) const = 0;
    //@}

        
};


#endif
