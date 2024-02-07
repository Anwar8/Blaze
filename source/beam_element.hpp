#ifndef BEAM_ELEMENT_HPP
#define BEAM_ELEMENT_HPP

#include <string>
#include <array>
#include <memory>
#include "maths_defaults.hpp"
#include "node.hpp"

/**
 * @brief basic cross-section containing geometry and material information and functions to get E, A, and I
 * 
 * @details This class currently only creates an elastic section with the equivalent properties of a 457 x 191 x 98
 * I-section made from steel. 
 */

// 457 x 191 x 98
class BasicSection {
    private:
        real E = 2.06e11; // Pa
        real A = 0.0125; // m^2
        real I = 0.0004570000; // m^4
    public:
        real get_E() {return E;}
        real get_A() {return A;}
        real get_I() {return I;}
};

/**
 * @brief shape function class that contains the local stiffness, freedoms, derivative, and order of DoFs
 * 
 * @details the shape function class is used to calculate the values of the shape function related matrices
 * for each row columns of the N and B. This requires knowing where along the element these values will be
 * calculated and how long the element is. Given section properties, this class will also calculate local
 * stiffness. Currently, this is done assuming an elastic non-changing local stiffness matrix.
 * 
 */
class BasicShapeFunction {
    private:
        mat k = make_xd_mat(6,6);
        mat N = make_xd_mat(2,6);
        mat B = make_xd_mat(2,6);
        std::vector<int> dof_map = {0, 2, 5};
    public:
        mat const get_k() const {return k;}
        mat const get_N() const {return N;}
        mat const get_B() const {return B;}
        std::vector<int> const get_dof_map() const {return dof_map;}
        /
        void calc_N(real x, real L);
        void calc_B(real x, real L);
        void calc_k(real L, BasicSection& sec);
};

/**
 * @brief 
 * 
 * creates transform matrices to shift and rotate elements
 * 
 * @details
 * 
 * calculates the local axis of a beam-column element, and then creates a transform
 * matrix T that will be used to transform the stiffness matrix.
 * 
 */

class BasicOrientation {
    private:
        coords local_x; /**< the x, y, z unit vectors for the beam element local x axis. */
        real length; /**< length of the beam-column element. */
        mat T = make_xd_mat(6,12); /**< The 6x12 T transformation matrix. */
        real alpha = 0.0; /**< angle between local and global coordinate systems. */
        real offset = 0.0; /**< offset of centroid along y axis. */
    public:
        /**
         * @brief sets offset and calculates length, local_x axis, alpha, and T matrix
         * 
         * @param nodes a size 2 array of shared ptr to the element nodes
         * @param sec_offset the offset value for the section
         * @param origin_x the global coord system x-axis unit vector
         */
        void evaluate(std::array<std::shared_ptr<Node>, 2> const & nodes, real sec_offset, coords const & origin_x)
        {
            offset = sec_offset;
            calc_length_local_x(nodes);
            calc_alpha(origin_x);
            calc_T();
        }
        /**
         * @brief calculates beam length and local_x components
         * 
         * @param nodes a size 2 array of shared ptr to the element nodes
         */
        void calc_length_local_x(std::array<std::shared_ptr<Node>, 2> const &  nodes) {
            local_x = (nodes[1]->get_coords() - nodes[0]->get_coords());
            length = local_x.norm();
            local_x /= length;
        }
        /**
         * @brief calcualte the angle between the local and global x axes
         * 
         * @details
         * The calculation is based on:
         * \f$ A\dot B = \abs{A}\abs{B}\cos{\theta}\f$
         * and then:
         * \f$ \theta = \arccos{\frac{A\dot B}{\abs{A}\abs{B}}} \f$
         * 
         * @param origin_x the global coord system x-axis unit vector
         */
        void calc_alpha(coords const& origin_x) {
           alpha = std::acos(origin_x.dot(local_x));
        }
        /**
         * @brief 
         * creates the T matrix
         * 
         * @details
         * the T matrix created is a 6x12 matrix with members only in the first
         * 6 rows. No idea why this matrix is not square.
         */
        void calc_T() {
            real c = std::cos(alpha);
            real s = std::sin(alpha);
            T(0,0) = c;
            T(0,2) = s;
            T(0,5) = offset;
            T(1,0) = -s;
            T(1,2) = c;
            T(2,5) = 1;
            T(3,6) = c;
            T(3,8) = s;
            T(3,11) = offset;
            T(4,6) = -s;
            T(4,8) = c;
            T(5,11) = 1;
        }
        mat get_T() {return T;}
};

/**
 * @brief a 2D beam element with 6 total freedoms: 1 rotation and two displacements
 * 
 * @details
 * This beam-column element has fewer freedoms than required in a 3D domain, and so the
 * transformation matrix T must map the element from 6 freedoms to the required 12 in 3D
 * domains.
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
        
        void move_nodes_up(real up) {
            for (auto node: nodes) {
                node->set_z(up);
            }
        }
        void set_d(vec new_disp) {local_d = new_disp;}

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
    // Copy assignment - i need to learn a bit more...
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
