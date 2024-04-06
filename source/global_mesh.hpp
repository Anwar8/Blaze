/**
 * @file global_mesh.hpp
 * @brief global mesh object, and functions related to creating model and operating on it
 * 
 */

#ifndef GLOBAL_MESH
#define GLOBAL_MESH

#include <memory>
#include <vector>
#include <tuple>
#include <string>
#include<Eigen/SparseLU>
#include<Eigen/SparseCholesky>
#include "gmsh.h"
#include "maths_defaults.hpp"
#include "node.hpp"
#include "beam_element.hpp"
#include "Izzuddin2DNonlinearBeam.hpp"

/**
 * @brief std vector of pairs each of which has an id and a 3-item coords vector.
 * 
 */
using gmsh_node_map = std::vector<std::pair<size_t, coords>>;

/**
 * @brief std vector of pairs each of which has an id and a vector of node ids related to the element id.
 * 
 */
using gmsh_elem_map = std::vector<std::pair<size_t, std::vector<size_t>>>;


/**
 * @brief contains the objects and functions operating on the global stiffness, and displacement and force vectors. Also contains the references to the nodes and elements.
 * 
 * @attention does not use the global matrices/vectors, so they should be removed as members.
 */
class GlobalMesh {
    private: 
        int nnodes = 0; /**< number of nodes in the mesh.*/
        int ndofs = 0; /**< number of DOFs in the mesh.*/
        int nelems = 0; /**< number of elements in the mesh.*/
        std::vector<std::shared_ptr<Node>> node_vector;  /**< a vector of shared ptrs referring to all the nodes in the problem.*/
        std::vector<std::shared_ptr<Basic2DBeamElement>> elem_vector; /**< a vector of shared ptrs referring to all the elements in the problem.*/
        // std::vector<std::shared_ptr<Izzuddin2DNonlinearBeam>> elem_vector; /**< a vector of shared ptrs referring to all the elements in the problem.*/

        spmat K;
        spvec P;
        vec U; 

    public:
        friend class Assembler;
        void open_mesh_file(std::string const mesh_file);
        /**
         * @brief uses gmsh API to read the nodes and populate a map with ids and coordinates.
         * 
         * @return gmsh_node_map a vector of pairs mapping each node ids and coordinates.
         */
        gmsh_node_map read_nodes();

        /**
         * @brief reads the elements from the mesh file and populates the map linking element id with its node ids.
         * 
         * @return gmsh_elem_map elem_map a vectpr of pairs mapping elem ids and corresponding node ids.
         */
        gmsh_elem_map read_elements();

        /**
         * @brief creates element objects following the \ref gmsh_elem_map object format and adds them to \ref elem_vector.
         * 
         * @param node_map a vector of pairs mapping each node ids and coordinates.
         */
        void make_nodes (gmsh_node_map node_map);

        /**
         * @brief creates element objects following the \ref gmsh_elem_map object format and adds them to \ref elem_vector.
         * 
         * @attention only one type of elements available now. Functionality limited to creating 2D beam-columns.
         * 
         * @param elem_map a vector of pairs mapping elem ids and corresponding node ids.
         */
        void make_elements (gmsh_elem_map elem_map);
        void close_mesh_file();

        /**
         * @brief populates the global_mesh members: \ref nnodes, \ref nelems, \ref node_vector, and \ref elem_vector.
         * 
         * @param mesh_file a string that is the mesh file name.
         */
        void setup_mesh(std::string const mesh_file);

        /**
         * @brief counts the active DOFs in the mesh by going over all the nodes and getting the number of active freedoms.
         * 
         */
        void count_dofs();
        /**
         * @brief calls the \ref Node::print_info function of all nodes and \ref Basic2DBeamElement::print_info of all elements.
         * 
         */
        void print_info();
        /**
         * @brief calls a node's \ref Node::fix_dof function.
         * 
         * @param id id of the node to which to apply the boundary condition.
         * @param dof DoF to apply the constraint.
         */
        void fix_node(int id, int dof);
        /**
         * @brief calls a node's \ref Node::add_nodal_load command and loads it.
         * 
         */
        void load_node(int id, int dof, real load);
        /**
         * @brief increments the nodal load of DoF dof of node id by load increment dP. Uses \ref Node::increment_nodal_load.
         * 
         * @param id 
         * @param dof 
         * @param dP
         */
        void increment_node_load(int id, int dof, real dP);

        /**
         * @brief keeps track of the history of a nodal DoF.
         * 
         * @param id 
         * @param dof 
         * @param history vector where nodal DoF history is tracked.
         */
        void track_nodal_dof(int id, int dof, std::vector<real>& history);

        /**
         * @brief loop over elements and call each of their \ref map_stiffness and \ref calc_K_global functions.
         * 
         */
        void calc_global_contributions() {
            std::cout << "Calc_global_contirbutions: There are " << ndofs << " active DoFs in the mesh." << std::endl;
            for (auto elem: elem_vector) 
            {   
                elem->map_stiffness();
                elem->calc_K_global();
            }
            for (auto node: node_vector)
            {
                node->compute_global_load_triplets();
            }
        }

        /**
         * @brief updates the state of each element after calculating global displacements.
         * 
         */
        void update_elements_states()
        {
            for (auto elem: elem_vector)
            {
                elem->update_state();
            }
        }
        /**
         * @brief prints the selected state of each element.
         * 
         */
        void print_elements_states(bool print_stresses = true, bool print_strains = false,
                                 bool print_nodal_disp = false, bool print_nodal_forces = false)
        {
            for (auto elem: elem_vector)
            {
                elem->print_element_state(print_stresses, print_strains, print_nodal_disp, print_nodal_forces);
            }
        }

        /**
         * @brief No idea why we have a template to solve for U here. 
         * @todo reconcile this function with the \ref BasicSolver object
         * 
         */
        void solve_for_U();
        int const get_num_elems() const {return nelems;}
};

#endif