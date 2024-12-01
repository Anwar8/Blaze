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
#include "ElementTypes.hpp"
#include "blaze_config.hpp"
#if INCLUDE_GMSH
    #include "gmsh.h"
#endif
#include "maths_defaults.hpp"
#include "node.hpp"
#include "BeamElementBaseClass.hpp"
#include "Linear2DBeamElement.hpp"
#include "Nonlinear2DBeamElement.hpp"
#include "Nonlinear2DPlasticBeamElement.hpp"
#include "basic_utilities.hpp"
#include "BeamColumnFiberSection.hpp"
#include "FrameMesh.hpp"

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
        std::vector<std::shared_ptr<ElementBaseClass>> elem_vector; /**< a vector of shared ptrs referring to all the elements in the problem.*/

        FrameMesh frame;
        
        // SectionBaseClass section; /**< a BasicSection object that is used by all elements in the mesh.*/
        std::unique_ptr<BasicSection> basic_section;
        std::unique_ptr<BeamColumnFiberSection> fiber_section;

        ElementType element_type;
        spmat K;
        spvec P;
        vec U; 

    public:
        friend class Assembler;
        #if INCLUDE_GMSH
        void open_mesh_file(std::string const mesh_file)
        {
            gmsh::initialize();
            gmsh::open(mesh_file);
        }
        /**
         * @brief uses gmsh API to read the nodes and populate a map with ids and coordinates.
         * 
         * @return NodeIdCoordsPairsVector a vector of pairs mapping each node ids and coordinates.
         */
        NodeIdCoordsPairsVector read_nodes()
        {
            std::vector<double> coord_vec;
            std::vector<double> parametricCoords;
            std::vector<std::size_t> nodeTags;

            gmsh::model::mesh::getNodes(nodeTags, coord_vec, parametricCoords);
            
            NodeIdCoordsPairsVector node_map;
            node_map.reserve(nodeTags.size());

            auto itr = coord_vec.begin();
            for (auto& tag : nodeTags)
            {
                node_map.push_back(std::make_pair(tag, coords(*itr, *(itr+1), *(itr+2))));
                itr += 3;
            }
            return node_map;    
        }

        /**
         * @brief reads the mesh file (gmsh format) and populates the node and element maps.
         * @param mesh_file a string that is the gmsh mesh file name and includes directory.
         * @param sect a \ref BasicSection object that is used to initialise the beam-column elements.
         * @return std::pair<NodeIdCoordsPairsVector, ElemIdNodeIdPairVector> a pair of node and element maps corresponding to a gmsh file.
         */
        std::pair<NodeIdCoordsPairsVector, ElemIdNodeIdPairVector> read_mesh_file(std::string const mesh_file, BasicSection sect) 
        {
            // section = sect;
            basic_section =  std::make_unique<BasicSection>(sect);
            open_mesh_file(mesh_file);
            NodeIdCoordsPairsVector node_map = read_nodes();
            ElemIdNodeIdPairVector elem_map = read_elements();
            close_mesh_file();
            return std::make_pair(node_map, elem_map);
        }
        
        /**
         * @brief reads the elements from the mesh file and populates the map linking element id with its node ids.
         * 
         * @return ElemIdNodeIdPairVector elem_map a vectpr of pairs mapping elem ids and corresponding node ids.
         */
        ElemIdNodeIdPairVector read_elements()
        {
            ElemIdNodeIdPairVector elem_map;
            std::vector<int> element_types;
            gmsh::model::mesh::getElementTypes(element_types);
            for (auto& elem_type: element_types)
            {
                std::string type_name;
                int num_nodes_per_elem;
                // throw-away variables
                int dim, order, num_primary_nodes;
                std::vector<double> local_nodal_coords;
                // get the information about the element type. in particular:
                // (1) number of nodes, and  (2) type name
                gmsh::model::mesh::getElementProperties(elem_type, type_name, dim, 
                                                        order, num_nodes_per_elem, local_nodal_coords, 
                                                        num_primary_nodes);
                
                // getting the elements and nodes for each element type
                // this is because each type has the same number of nodes
                std::cout << "Mapping element type: " << type_name << std::endl;
                std::vector<std::size_t> elem_tags, node_tags;
                gmsh::model::mesh::getElementsByType(elem_type, elem_tags, node_tags);

                std::vector<std::size_t> element_nodes;
                element_nodes.reserve(num_nodes_per_elem);
                auto node_itr = node_tags.begin() ;
                for (auto& elem_tag: elem_tags)
                {
                    element_nodes.clear();
                    for (int i = 0; i < num_nodes_per_elem; ++i)
                    {
                        element_nodes.push_back(*(node_itr + i));
                        node_itr += i;
                    }
                    ++node_itr;
                    elem_map.push_back(std::make_pair(elem_tag, element_nodes));
                    std::cout << "added element: " << elem_tag << " with nodes: ";
                    print_container<std::vector<std::size_t>>(element_nodes);
                }
            }
            return elem_map;

            void GlobalMesh::close_mesh_file()
            {
            gmsh::finalize();
            }
        }
        #endif

        void create_frame_mesh(int nbays, int nfloors, real bay_length, real floor_height, int beam_divisions, int column_divisions, ElementType elem_type, BeamColumnFiberSection& sect)
        {
            frame = FrameMesh(nbays, nfloors, bay_length, floor_height, beam_divisions, column_divisions);
            element_type = elem_type;
            fiber_section = std::make_unique<BeamColumnFiberSection>(sect);

            NodeIdCoordsPairsVector node_map = frame.get_node_coords_pairs();
            ElemIdNodeIdPairVector elem_map = frame.map_elements_to_nodes();
            std::cout << "-----------------------------------------------------" << std::endl;
            read_node_map(node_map);
            std::cout << "-----------------------------------------------------" << std::endl;
            read_element_map(elem_map);
            std::cout << "-----------------------------------------------------" << std::endl;
            setup_mesh(node_map, elem_map);
        }

        FrameMesh get_frame() {return frame;}

        /**
         * @brief creates a line mesh map with a given number of divisions and end coordinates of the line. Does NOT take a gmsh mesh file.
         * 
         * @tparam CoordsContainer a container that has the coordinates of the end points of the line. Needs compatible with STL iterators.
         * @param divisions number of divisions to break the line into.
         * @param pts_coords the coordinates of the end points of the line.
         * @param elem_type an enum referring to the type of element that the mesh will include.
         * @param sect a \ref BasicSection object that is used to initialise the beam-column elements.
         * @return std::pair<NodeIdCoordsPairsVector, ElemIdNodeIdPairVector> the node_map and elem_map of the line mesh. 
         * @warning assumes mapping takes place from node and element ids = 1. There is no checking for conflicting ids, and nothing to reduce bandwidth!
         */
        template <typename CoordsContainer>
        std::pair<NodeIdCoordsPairsVector, ElemIdNodeIdPairVector> map_a_line_mesh(unsigned divisions, CoordsContainer pts_coords)
        {
            if (pts_coords.size() != 2)
            {
                std::cout << "Error: pts_coords must have 2 elements." << std::endl;
                exit(1);
            }

            NodeIdCoordsPairsVector node_map;
            ElemIdNodeIdPairVector elem_map;
            node_map.reserve(divisions + 1);
            elem_map.reserve(divisions);


            coords delta_xyz = (pts_coords[1] - pts_coords[0])/divisions;

            for (size_t i = 0; i < divisions + 1; ++i)
            {
                node_map.push_back(std::make_pair(i + 1, pts_coords[0] + i*delta_xyz));
            }

            for (size_t i = 0; i < divisions; ++i)
            {
                elem_map.push_back(std::make_pair(i + 1, std::vector<size_t>{i + 1, i + 2}));
            }     
            return std::make_pair(node_map, elem_map);
        }
        /**
         * @brief creates a line mesh map with a given number of divisions and end coordinates of the line. Does NOT take a gmsh mesh file. 
         * @param elem_type an enum referring to the type of element that the mesh will include.
         * @param sect a \ref BeamColumnFiberSection object that is used to initialise the beam-column elements.
         * 
        **/
        void create_line_mesh(int divisions, std::vector<coords> end_coords, ElementType elem_type, BeamColumnFiberSection& sect)
        {
            element_type = elem_type;
            fiber_section = std::make_unique<BeamColumnFiberSection>(sect);
            std::pair<NodeIdCoordsPairsVector, ElemIdNodeIdPairVector> mesh_maps = map_a_line_mesh(divisions, end_coords);
            setup_mesh(mesh_maps.first, mesh_maps.second);
        }
        /**
         * @brief creates a line mesh map with a given number of divisions and end coordinates of the line. Does NOT take a gmsh mesh file. 
         * @param elem_type an enum referring to the type of element that the mesh will include.
         * @param sect a \ref BasicSection object that is used to initialise the beam-column elements.
         * 
        **/
        void create_line_mesh(int divisions, std::vector<coords> end_coords, ElementType elem_type, BasicSection& sect)
        {
            element_type = elem_type;
            basic_section = std::make_unique<BasicSection>(sect);
            std::pair<NodeIdCoordsPairsVector, ElemIdNodeIdPairVector> mesh_maps = map_a_line_mesh(divisions, end_coords);
            setup_mesh(mesh_maps.first, mesh_maps.second);
        }


        /**
         * @brief creates element objects following the \ref ElemIdNodeIdPairVector object format and adds them to \ref elem_vector.
         * 
         * @param node_map a vector of pairs mapping each node ids and coordinates.
         */
        void make_nodes (NodeIdCoordsPairsVector node_map)
        {
            for (auto& node_data : node_map)
            {
                node_vector.push_back(std::make_shared<Node>(node_data.first, node_data.second));
            }
        }

        /**
         * @brief creates element objects following the \ref ElemIdNodeIdPairVector object format and adds them to \ref elem_vector.
         * 
         * @attention only one type of elements available now. Functionality limited to creating 2D beam-columns.
         * 
         * @param elem_map a vector of pairs mapping elem ids and corresponding node ids.
         */
        void make_elements (ElemIdNodeIdPairVector elem_map)
        {
            std::vector<std::shared_ptr<Node>> elem_nodes;
            elem_nodes.reserve(2);
            for (auto& element_data : elem_map)
            {
                elem_nodes.clear();
                for (auto& node_id: element_data.second)
                {
                    auto node = get_id_iterator<std::vector<std::shared_ptr<Node>>::iterator, std::vector<std::shared_ptr<Node>>>(node_id, node_vector);
                    elem_nodes.push_back(*node);
                }

                switch (element_type)
                {
                case LinearElastic:
                    elem_vector.emplace_back(std::make_shared<Linear2DBeamElement>(element_data.first, elem_nodes, *basic_section));
                    break;
                case NonlinearElastic:
                    elem_vector.emplace_back(std::make_shared<Nonlinear2DBeamElement>(element_data.first, elem_nodes, *basic_section));
                    break;
                case NonlinearPlastic:
                    elem_vector.emplace_back(std::make_shared<Nonlinear2DPlasticBeamElement>(element_data.first, elem_nodes, *fiber_section));
                    break;
                default:
                    std::cout << "Incorrect element_type: " << element_type << "; should be 0. LinearElastic, 1. NonlinearElastic, or 2. NonlinearPlastic." << std::endl;
                    exit(1);
                }
            }       
        }
        void close_mesh_file();


        
        /**
         * @brief populates the global_mesh members: \ref nnodes, \ref nelems, \ref node_vector, and \ref elem_vector based on mesh (node and element) maps.
         * 
         * @param node_map a vector of pairs mapping each node ids and coordinates.
         * @param elem_map a vector of pairs mapping elem ids and corresponding 2 node ids.
         */
        void setup_mesh(NodeIdCoordsPairsVector node_map, ElemIdNodeIdPairVector elem_map)
        {
            nnodes = node_map.size();
            nelems = elem_map.size();
            node_vector.clear();
            node_vector.reserve(nnodes);
            elem_vector.clear();
            elem_vector.reserve(nelems);
            make_nodes(node_map);
            make_elements(elem_map);
            std::sort(node_vector.begin(), node_vector.end());
            std::sort(elem_vector.begin(), elem_vector.end());
            count_dofs();            
        }

        /**
         * @brief counts the active DOFs in the mesh by going over all the nodes and getting the number of active freedoms.
         * 
         */
        void count_dofs()
        {
            ndofs = 0;
            for (auto& node: node_vector)
            {
                node->set_nz_i(ndofs);
                if (VERBOSE)
                {
                    std::cout << "Node " << node->get_id() << " has nz_i = " << ndofs << std::endl;
                }
                ndofs += node->get_ndof();
            }
        }
        /**
         * @brief calls the \ref Node::print_info function of all nodes and \ref Basic2DBeamElement::print_info of all elements.
         * 
         */
        void print_info()
        {
            std::cout << "Mesh contains " << nelems << " elements and " << nnodes << " nodes." << std::endl;
            for (auto& node: node_vector)
            {
                node->print_info();
            }
            for (auto& elem: elem_vector)
            {
                elem->print_info();
            }            
        }
        /**
         * @brief calls a node's \ref Node::fix_dof function.
         * 
         * @param id id of the node to which to apply the boundary condition.
         * @param dof DoF to apply the constraint.
         */
        void fix_node(int id, int dof)
        {
            auto node_it = get_id_iterator<std::vector<std::shared_ptr<Node>>::iterator, std::vector<std::shared_ptr<Node>>>(id, node_vector);
            if (dof < 0)
            {
                std::cout << "Fixing all DoFs of node " << id << std::endl;
                (*node_it)->fix_all_dofs();
                (*node_it)->print_inactive_dofs();
            } else {
                (*node_it)->fix_dof(dof);
            }    
        }
        /**
         * @brief calls a node's \ref Node::add_nodal_load command and loads it.
         * 
         */
        void load_node(int id, int dof, real load)
        {
            get_node_by_id(id)->add_nodal_load(load, dof);
        }
        /**
         * @brief increments the nodal load of DoF dof of node id by load increment dP. Uses \ref Node::increment_nodal_load.
         * 
         * @param id 
         * @param dof 
         * @param dP
         */
        void increment_node_load(int id, int dof, real dP)
        {
            get_node_by_id(id)->increment_nodal_load(dP, dof);
        }

        /**
         * @brief keeps track of the history of a nodal DoF.
         * 
         * @param id 
         * @param dof 
         * @param history vector where nodal DoF history is tracked.
         */
        void track_nodal_dof(int id, int dof, std::vector<real>& history)
        {
            auto node_it = get_id_iterator<std::vector<std::shared_ptr<Node>>::iterator, std::vector<std::shared_ptr<Node>>>(id, node_vector);
            std::array<real, 6> nodal_displacements = (*node_it)->get_nodal_displacements();
            history.push_back(nodal_displacements[dof]);
        }

        /**
         * @brief get node shared_ptr by id.
         * @param id id of the node to get from the \ref node_vector.
         * @return std::shared_ptr<Node> a shared pointer to the node with the given id.
         */
        std::shared_ptr<Node> get_node_by_id(int id)
        {
            auto node_it = get_id_iterator<std::vector<std::shared_ptr<Node>>::iterator, std::vector<std::shared_ptr<Node>>>(id, node_vector);
            return *node_it;
        }

        /**
         * @brief maps the stiffness of each element in the \ref elem_vector by calling \ref BeamElementCommonInterface::map_stiffness.
         * @details this function should be called after defining the constraints on the nodes, but before actually beginning th solution procedure.
         */
        void map_element_stiffnesses()
        {
            for (auto& elem: elem_vector)
            {
                elem->map_stiffness();
            }
        }

        /**
         * @brief loop over elements and call each of their \ref map_stiffness and \ref calc_global_stiffness_triplets functions.
         * 
         */
        void calc_global_contributions() {
            if (VERBOSE)
            {
                std::cout << "Calc_global_contirbutions: There are " << ndofs << " active DoFs in the mesh." << std::endl;
            }
            for (auto& elem: elem_vector) 
            {   
                // elem->update_state(); // Should not update state if we are calling update_state explicitly on its own!
                // elem->map_stiffness();
                elem->calc_global_stiffness_triplets();
            }
            for (auto& node: node_vector)
            {
                if (VERBOSE)
                {
                    std::cout << "Computing global load triplets for node " << node->get_id() << std::endl;
                }
                node->compute_global_load_triplets();
            }
        }

        /**
         * @brief checks the nodal loads of all nodes in the mesh using the \ref Node::check_loads function, which prints an error if a DoF is both loaded and restrained.
         * 
         */
        void check_nodal_loads()
        {
            for (auto& node: node_vector)
            {
                node->check_loads();
            }
        }
        /**
         * @brief updates the state of each element after calculating global displacements.
         * 
         */
        void update_elements_states()
        {
            #pragma omp parallel for
            for (auto& elem: elem_vector)
            {
                elem->update_state();
            }
        }

        void update_element_sections_starting_states()
        {
            for (auto& elem: elem_vector)
            {
                elem->update_section_starting_state();
            }
        }
        /**
         * @brief prints the selected state of each element.
         * 
         */
        void print_elements_states(bool print_nodal_disp = false, bool print_strains = false, 
                                            bool print_stresses = true, bool print_nodal_forces = false)
        {
            for (auto& elem: elem_vector)
            {
                elem->print_element_state(print_nodal_disp, print_strains, print_stresses, print_nodal_forces);
            }
        }

        /**
         * @brief No idea why we have a template to solve for U here. 
         * @todo reconcile this function with the \ref BasicSolver object
         * 
         */
        void solve_for_U()
        {
            Eigen::SparseLU<spmat> solver;
            // Compute the ordering permutation vector from the structural pattern of A
            solver.analyzePattern(K); 
            // Compute the numerical factorization 
            solver.factorize(K); 
            //Use the factors to solve the linear system 
            
            if (solver.info() == Eigen::Success)
            {
                std::cout << "Factorisation successfull." << std::endl;
            } else {
                std::cout << "ERROR: Factorisation unsuccessfull! Matrix is:" << std::endl;
                // convert to dense matrix to print correctly
                std::cout << Eigen::MatrixXd(K) << std::endl;       
                std::exit(1);
            }
            U = solver.solve(P); 
            std::cout << "The solution is:" << std::endl << U << std::endl;    
        }
        int const get_num_elems() const {return nelems;}
        int const get_num_nodes() const {return nnodes;}
};

#endif
