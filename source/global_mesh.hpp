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
#include "main.hpp"
#include "gmsh.h"
#include "maths_defaults.hpp"
#include "node.hpp"
#include "beam_element.hpp"
#include "Izzuddin2DNonlinearBeam.hpp"
#include "BeamElementBaseClass.hpp"

/**
 * @brief std vector of pairs each of which has an id and a 3-item coords vector.
 * 
 */
using NodeIdCoordsPairsVector = std::vector<std::pair<size_t, coords>>;

/**
 * @brief std vector of pairs each of which has an id and a vector of node ids related to the element id.
 * 
 */
using ElemIdNodeIdPairVector = std::vector<std::pair<size_t, std::vector<size_t>>>;


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
        #if (ELEM == 1 || ELEM == 2) 
            std::vector<std::shared_ptr<Basic2DBeamElement>> elem_vector; /**< a vector of shared ptrs referring to all the elements in the problem.*/
        #elif (ELEM == 3)
            std::vector<std::shared_ptr<BeamElementBaseClass>> elem_vector; /**< a vector of shared ptrs referring to all the elements in the problem.*/
        #else 
            std::cout << "Incorrect ELEM: " << ELEM << "; should be 1. OLD, 2. IZDN, or 3. LBE." << std::endl;
            exit(1);
        #endif
        

        spmat K;
        spvec P;
        vec U; 

    public:
        friend class Assembler;
        void open_mesh_file(std::string const mesh_file);
        /**
         * @brief uses gmsh API to read the nodes and populate a map with ids and coordinates.
         * 
         * @return NodeIdCoordsPairsVector a vector of pairs mapping each node ids and coordinates.
         */
        NodeIdCoordsPairsVector read_nodes();


        /**
         * @brief reads the elements from the mesh file and populates the map linking element id with its node ids.
         * 
         * @return ElemIdNodeIdPairVector elem_map a vectpr of pairs mapping elem ids and corresponding node ids.
         */
        ElemIdNodeIdPairVector read_elements();

        /**
         * @brief creates a line mesh map with a given number of divisions and end coordinates of the line. Does NOT take a gmsh mesh file.
         * 
         * @tparam CoordsContainer a container that has the coordinates of the end points of the line. Needs compatible with STL iterators.
         * @param divisions number of divisions to break the line into.
         * @param end_coords the coordinates of the end points of the line.
         * @return std::pair<NodeIdCoordsPairsVector, ElemIdNodeIdPairVector> the node_map and elem_map of the line mesh. 
         * @warning assumes mapping takes place from node and element ids = 1. There is no checking for conflicting ids, and nothing to reduce bandwidth!
         */
        template <typename CoordsContainer>
        std::pair<NodeIdCoordsPairsVector, ElemIdNodeIdPairVector> map_a_line_mesh(unsigned divisions, CoordsContainer end_coords)
        {
            if (end_coords.size() != 2)
            {
                std::cout << "Error: end_coords must have 2 elements." << std::endl;
                exit(1);
            }

            NodeIdCoordsPairsVector node_map;
            ElemIdNodeIdPairVector elem_map;
            node_map.reserve(divisions + 1);
            elem_map.reserve(divisions);


            coords delta_xyz = (end_coords[1] - end_coords[0])/divisions;

            for (size_t i = 0; i < divisions + 1; ++i)
            {
                node_map.push_back(std::make_pair(i + 1, end_coords[0] + i*delta_xyz));
            }

            for (size_t i = 0; i < divisions; ++i)
            {
                elem_map.push_back(std::make_pair(i + 1, std::vector<size_t>{i + 1, i + 2}));
            }     
            return std::make_pair(node_map, elem_map);
        }

        /**
         * @brief creates element objects following the \ref ElemIdNodeIdPairVector object format and adds them to \ref elem_vector.
         * 
         * @param node_map a vector of pairs mapping each node ids and coordinates.
         */
        void make_nodes (NodeIdCoordsPairsVector node_map);

        /**
         * @brief creates element objects following the \ref ElemIdNodeIdPairVector object format and adds them to \ref elem_vector.
         * 
         * @attention only one type of elements available now. Functionality limited to creating 2D beam-columns.
         * 
         * @param elem_map a vector of pairs mapping elem ids and corresponding node ids.
         */
        void make_elements (ElemIdNodeIdPairVector elem_map);
        void close_mesh_file();

        /**
         * @brief reads the mesh file (gmsh format) and populates the node and element maps.
         * 
         * @param mesh_file a string that is the gmsh mesh file name and includes directory.
         * @return std::pair<NodeIdCoordsPairsVector, ElemIdNodeIdPairVector> a pair of node and element maps corresponding to a gmsh file.
         */
        std::pair<NodeIdCoordsPairsVector, ElemIdNodeIdPairVector> read_mesh_file(std::string const mesh_file) 
        {
            open_mesh_file(mesh_file);
            NodeIdCoordsPairsVector node_map = read_nodes();
            ElemIdNodeIdPairVector elem_map = read_elements();
            close_mesh_file();
            return std::make_pair(node_map, elem_map);
        }
        
        /**
         * @brief populates the global_mesh members: \ref nnodes, \ref nelems, \ref node_vector, and \ref elem_vector based on mesh (node and element) maps.
         * 
         * @param node_map a vector of pairs mapping each node ids and coordinates.
         * @param elem_map a vector of pairs mapping elem ids and corresponding 2 node ids.
         */
        void setup_mesh(NodeIdCoordsPairsVector node_map, ElemIdNodeIdPairVector elem_map);

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
         * @brief loop over elements and call each of their \ref map_stiffness and \ref calc_global_stiffness_triplets functions.
         * 
         */
        void calc_global_contributions() {
            std::cout << "Calc_global_contirbutions: There are " << ndofs << " active DoFs in the mesh." << std::endl;
            for (auto elem: elem_vector) 
            {   
                elem->update_state();
                elem->map_stiffness();
                elem->calc_global_stiffness_triplets();
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
        void print_elements_states(bool print_nodal_disp = false, bool print_strains = false, 
                                            bool print_stresses = true, bool print_nodal_forces = false)
        {
            for (auto elem: elem_vector)
            {
                elem->print_element_state(print_nodal_disp, print_strains, print_stresses, print_nodal_forces);
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