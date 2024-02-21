/**
 * @file beam_element.hpp
 * @brief definitions for basic beam-column elements.
 * 
 */

#ifndef BEAM_ELEMENT_HPP
#define BEAM_ELEMENT_HPP

#include <string>
#include <array>
#include <memory>
#include "maths_defaults.hpp"
#include "basic_section.hpp"
#include "basic_shape_function.hpp"
#include "basic_orientation.hpp"
#include "node.hpp"

/**
 * @brief a 2D beam element with 6 total freedoms: 1 rotation and two displacements
 * 
 * @details
 * This beam-column element has fewer freedoms than required in a 3D domain, and so the
 * transformation matrix T must map the element from 6 freedoms to the required 12 in 3D
 * domains.
 * 
 * @todo 1. add functionality to apply distributed mechanical loading
 *       2. add functionality to apply thermal gradient and uniform temperature
 *       3. add functionality to apply non-uniform thermal loading interpolation (LOW PRIORITY)
 * 
 */
class Basic2DBeamElement {
    private:
        unsigned id = 0; /**< unique identifier for the element.*/
        std::string const elem_type = "beam-column"; /**< string that represents the type of the element.*/
        int const ndofs = 3; /**< number of freedoms at each node.*/
        int const nnodes = 2; /**< number of nodes.*/
        std::array<std::shared_ptr<Node>, 2> nodes; /**< an array of size 2 that holds the shared ptrs to the nodes.*/
        
        std::vector<int> global_dof_map; /**< FORGOT - appears to have been deprecated.*/


        /**
         * @brief maps local stiffness contributions to their global positions in the stiffness matrix
         * 
         * @details a vector of std array of size 4. the first two indices of the array refer to the
         * transformed local stiffness matrix indices, and the last two refer to the indices where that
         * local stiffness would go in the global stiffness matrix. The size of the std vector is dependent on the
         * number of nodes and which DOFs are active/not fixed. See \ref map_stiffness for details.
         */
        std::vector<std::array<int,4>> stiffness_map; 

        BasicSection section; /**< the section for the beam-column element.*/
        BasicShapeFunction shape_func; /**< the shape function for the beam-column element.*/
        BasicOrientation orient; /**< the orientation object for the beam-column element.*/
        real length = 0.0; /**< the length for the beam-column element - to be calculated by the orientation object.*/

        vec local_d = make_xd_vec(6); /**< local displacements for all freedoms.*/
        vec local_f = make_xd_vec(6); /**< local forces corresponding to all freedoms.*/
        vec local_eps = make_xd_vec(2); /**< local strains corresponding to all freedoms.*/

        std::vector<spnz> K_global; /**< the global contributions of the element to the global stiffness - made as sparse matrix contributions that would be gatehred to create the global sparse matrix*/

    public:
        /**
         * @brief Construct a new Basic 2D Beam Element object (Default)
         * 
         * @details Kept, most likely, for copying...otherwise, I don't see how this is
         * going to work. Probably need to create a copy constructor as there are some
         * fishy pointer stuff going on
         * 
         * @todo deprecate this constructor
         */
        Basic2DBeamElement();

        /**
         * @brief Construct a new Basic 2D Beam Element object when given two nodes, and calculate its length
         * 
         * @details
         * Definitely needs to be deprecated as it simply passes the default element id
         * which is 0 (or unintialized) to the nodes
         * 
         * @todo deprecate this constructor
         * 
         * @param node_1 a shared ptr referencing the first node
         * @param node_2 a shared ptr referencing the second node
         */
        Basic2DBeamElement(std::shared_ptr<Node>& node_1, std::shared_ptr<Node>& node_2);

        /**
         * @brief Construct a new Basic 2D Beam Element object when given two nodes and an id, and calculate its length
         * 
         * @details 
         * Nodes are given the ids of the elements connected to them to facilitate node deletion
         * in future. Right now, this is just creating extra mess.
         * 
         * @param id unique identifier for the element; will be passed to the nodes
         * @param node_1 a shared ptr referencing the first node
         * @param node_2 a shared ptr referencing the second node
         */
        Basic2DBeamElement(int id, std::shared_ptr<Node>& node_1, std::shared_ptr<Node>& node_2);
        /**
         * @brief Construct a new Basic 2D Beam Element object from a vector of nodes and an id, and calculate its length
         * 
         * @details 
         * Actually easier to use as we don't have to give it two nodes and instead can easily give
         * a container of nodes. Will check number of nodes is \ref nnodes to match the element requirements.
         * 
         * @todo    deprecate in favour of template that accepts any container type
         * 
         * @param given_id unique identifier for the element; will be passed to the nodes
         * @param in_nodes a std vector of shared pointers to node objects
         */
        Basic2DBeamElement(int given_id, std::vector<std::shared_ptr<Node>>& in_nodes);
        
        /**
         * @brief Construct a new Basic 2D Beam Element object and calculate its length
         * 
         * @tparam Container any type of std container that has a std::size and built-in iterators
         * @param given_id unique identifier for the element; will be passed to the nodes
         * @param in_nodes a container of shared pointers to node objects
         */
        template<typename Container>
        Basic2DBeamElement(int given_id, Container& in_nodes) {
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
        }
        void print_info();
        void calc_length();
        void calc_N(real x);
        void calc_B(real x);
        void calc_k();
        void calc_T(real sec_offset = 0.0, coords origin_x = {1.0, 0.0, 0.0});
        void calc_eps();

        /**
         * @brief calculates the global stiffness contribution of the local element and populates K_global
         * 
         * @details first, the freedoms are mapped to the right size by pre- and post-multiplying by the T matrix
         * After that, \ref stiffness_map is used to map where these contributions would go in the global stiffness
         * matrix. So, this function will populate \ref K_global with sparse matrix notation
         * 
         */
        void calc_K_global();

        /**
         * @brief populates \ref stiffness_map considering active and inactive DOFs for each node of the element
         * 
         * @details see function \ref calc_K_global, and variables \ref stiffness_map, and \ref K_global. 
         * 
         * @todo REALLY needs to be revisited. attempt to rewrite this function so it does the following:
         *  1. gets all the contribution without worrying about active or not
         *  2. if a contribution is inactive then that contribution is zeroed AND
         *  3. zeroed contributions are not added to \ref K_global
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
        std::vector<int> get_global_dof_map() {return global_dof_map;}
        std::vector<spnz> get_K_global() {return K_global;}

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
    /**
     * @brief Copy assignment - i need to learn a bit more to make sure this actually works.
     * 
     * @param other 
     * @return Basic2DBeamElement& 
     */
    Basic2DBeamElement(const Basic2DBeamElement& other) {
        id = other.id;
        nodes = other.nodes;
        section = other.section;
        shape_func = other.shape_func;
        length = other.length;
        local_d = other.local_d;
        local_f = other.local_f;
        local_eps = other.local_eps;
    }
    /**
     * @brief Copy assignment - i need to learn a bit more...
     * 
     * @param other 
     * @return Basic2DBeamElement& 
     */
    Basic2DBeamElement& operator=(const Basic2DBeamElement& other) {
        if (this != &other) {
            Basic2DBeamElement temp(other); // use copy constructor to create temporary object

            std::swap(id, temp.id);
            std::swap(nodes, temp.nodes);
            std::swap(section, temp.section);
            std::swap(shape_func, temp.shape_func);
            std::swap(length, temp.length);
            std::swap(local_d, temp.local_d);
            std::swap(local_f, temp.local_f);
            std::swap(local_eps, temp.local_eps);
        }
        return *this;
    }
};


#endif
