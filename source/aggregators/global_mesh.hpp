/**
 * @file global_mesh.hpp
 * @brief global mesh object, and functions related to creating model and operating on it
 * 
 */

#ifndef GLOBAL_MESH
#define GLOBAL_MESH

// debugging!! Remove!
#include <unistd.h>


#include <memory>
#include <vector>
#include <tuple>
#include <string>
#include <map>
#include <set>
#include <algorithm>
#include <Eigen/SparseLU>
#include <Eigen/SparseCholesky>
#ifdef WITH_MPI
    #include <mpi.h>
    #include <tpetra_wrappers.hpp>
#endif
#ifdef KOKKOS
    #include <Kokkos_Core.hpp>
#endif
#include "ElementTypes.hpp"
#include "blaze_config.hpp"
#if INCLUDE_GMSH
    #include "gmsh.h"
#endif
#include "maths_defaults.hpp"
#include "MPIWrappers.hpp"
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
using NodeIdCoordsPairsVector = std::vector<std::pair<unsigned, coords>>;

/**
 * @brief std vector of pairs each of which has an element id and a vector of node ids related to the element id.
 * 
 */
using ElemIdNodeIdPairVector = std::vector<std::pair<unsigned, std::vector<unsigned>>>;


/**
 * @brief contains the objects and functions operating on the global stiffness, and displacement and force vectors. Also contains the references to the nodes and elements.
 * 
 * @attention does not use the global matrices/vectors, so they should be removed as members.
 */
class GlobalMesh {
    private: 
        int rank = 0; /**< the rank on which this GlobalMesh object is created.*/
        int num_ranks = 1; /**< the total number of ranks running on the MPI processes. */
        int nnodes = 0; /**< number of nodes in the mesh.*/
        int ndofs = 0; /**< number of DOFs in the mesh.*/
        int nelems = 0; /**< number of elements in the mesh.*/
        int rank_nnodes = 0; /**< number of nodes on current rank.*/
        int rank_interface_nnodes = 0; /**< number of interface nodes on current rank.*/
        int rank_ndofs = 0; /**< number of DOFs on current rank.*/
        int rank_nelems = 0; /**< number of elements on current rank.*/
        unsigned rank_starting_node_id = 1; /**< number at which the id of the first node on this rank starts.*/
        int rank_starting_nz_i = 0; /**< number at which the DoF count starts on this rank. */
        int max_num_stiffness_contributions = 0; /**< finds the maximum number of contributions made to a row of the stiffness matrix */
        
        std::vector<int> ranks_ndofs; /**< a vector that contains the \ref ndofs of each rank. Needs to be communicated with MPI across ranks.*/
        std::vector<unsigned> ranks_nnodes; /**< a vector that contains the \ref rank_nnodes of each rank. Needs to be communicated with MPI across ranks.*/
        std::vector<std::shared_ptr<Node>> node_vector;  /**< a vector of shared ptrs referring to all the nodes on the current rank.*/
        std::vector<std::shared_ptr<Node>> interface_node_vector;  /**< a vector of shared ptrs referring to the interface nodes in the problem.*/
        std::vector<std::shared_ptr<ElementBaseClass>> elem_vector; /**< a vector of shared ptrs referring to all the elements in the problem.*/

        std::set<unsigned> interface_node_id_set_on_rank; //**<a std::set of node ids for those that are interfaces on this rank. */
        std::set<unsigned> node_id_set_owned_by_rank; //**<a std::set of node ids for those that are owned by this rank. */


        std::map<int, std::set<unsigned>> wanted_by_neighbour_rank_node_id_map;/**<a std::map that maps neighbour rank to a set of node IDs that this neighbour wants from the current rank. */
        std::map<int, std::set<unsigned>> wanted_from_neighbour_rank_node_id_map; /**<a std::map that maps neighbour rank to a set of node IDs that the current rank wants from this neighbour. */


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
            std::vector<std::unsigned> nodeTags;

            gmsh::model::mesh::getNodes(nodeTags, coord_vec, parametricCoords);
            
            NodeIdCoordsPairsVector nodes_coords_vector;
            nodes_coords_vector.reserve(nodeTags.size());

            auto itr = coord_vec.begin();
            for (auto& tag : nodeTags)
            {
                nodes_coords_vector.push_back(std::make_pair(tag, coords(*itr, *(itr+1), *(itr+2))));
                itr += 3;
            }
            return nodes_coords_vector;    
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
            NodeIdCoordsPairsVector nodes_coords_vector = read_nodes();
            ElemIdNodeIdPairVector elem_nodes_vector = read_elements();
            close_mesh_file();
            return std::make_pair(nodes_coords_vector, elem_nodes_vector);
        }
        
        /**
         * @brief reads the elements from the mesh file and populates the map linking element id with its node ids.
         * 
         * @return ElemIdNodeIdPairVector elem_nodes_vector a vectpr of pairs mapping elem ids and corresponding node ids.
         */
        ElemIdNodeIdPairVector read_elements()
        {
            ElemIdNodeIdPairVector elem_nodes_vector;
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
                std::vector<std::unsigned> elem_tags, node_tags;
                gmsh::model::mesh::getElementsByType(elem_type, elem_tags, node_tags);

                std::vector<std::unsigned> element_nodes;
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
                    elem_nodes_vector.push_back(std::make_pair(elem_tag, element_nodes));
                    std::cout << "added element: " << elem_tag << " with nodes: ";
                    print_container<std::vector<std::unsigned>>(element_nodes);
                }
            }
            return elem_nodes_vector;

            void GlobalMesh::close_mesh_file()
            {
            gmsh::finalize();
            }
        }
        #endif

        /**
         * @brief sets the value for \ref rank and \ref num_ranks based on MPI environment.
         */
        void initialise_mpi_variables()
        {
            get_my_rank(rank);
            get_num_ranks(num_ranks);
        }

        int get_mesh_num_ranks() {return num_ranks;}
        int get_mesh_rank() {return rank;}


        void create_frame_mesh(int nbays, int nfloors, real bay_length, real floor_height, int beam_divisions, int column_divisions, ElementType elem_type, BeamColumnFiberSection& sect)
        {
            frame = FrameMesh(nbays, nfloors, bay_length, floor_height, beam_divisions, column_divisions);
            element_type = elem_type;
            fiber_section = std::make_unique<BeamColumnFiberSection>(sect);

            NodeIdCoordsPairsVector nodes_coords_vector = frame.get_node_coords_pairs();
            ElemIdNodeIdPairVector elem_nodes_vector = frame.map_elements_to_nodes();
            if (VERBOSE)
            {
                std::cout << "-----------------------------------------------------" << std::endl;
                read_nodes_coords_vector(nodes_coords_vector);
                std::cout << "-----------------------------------------------------" << std::endl;
                read_element_map(elem_nodes_vector);
                std::cout << "-----------------------------------------------------" << std::endl;
            }
            setup_mesh(nodes_coords_vector, elem_nodes_vector);
        }

        void create_distributed_frame_mesh(int nbays, int nfloors, real bay_length, real floor_height, int beam_divisions, int column_divisions, ElementType elem_type, BeamColumnFiberSection& sect)
        {
            initialise_mpi_variables();
            frame = FrameMesh(nbays, nfloors, bay_length, floor_height, beam_divisions, column_divisions);
            element_type = elem_type;
            fiber_section = std::make_unique<BeamColumnFiberSection>(sect);

            NodeIdCoordsPairsVector nodes_coords_vector = frame.get_node_coords_pairs();
            ElemIdNodeIdPairVector elem_nodes_vector = frame.map_elements_to_nodes();
            if (VERBOSE)
            {
                std::cout << "-----------------------------------------------------" << std::endl;
                read_nodes_coords_vector(nodes_coords_vector);
                std::cout << "-----------------------------------------------------" << std::endl;
                read_element_map(elem_nodes_vector);
                std::cout << "-----------------------------------------------------" << std::endl;
            }
            setup_distributed_mesh(nodes_coords_vector, elem_nodes_vector);
        }

        void create_distributed_frame_mesh(int nbays, int nfloors, real bay_length, real floor_height, int beam_divisions, int column_divisions, ElementType elem_type, BasicSection& sect)
        {
            initialise_mpi_variables();
            frame = FrameMesh(nbays, nfloors, bay_length, floor_height, beam_divisions, column_divisions);
            element_type = elem_type;
            set_basic_section(sect);

            NodeIdCoordsPairsVector nodes_coords_vector = frame.get_node_coords_pairs();
            ElemIdNodeIdPairVector elem_nodes_vector = frame.map_elements_to_nodes();
            if (VERBOSE)
            {
                std::cout << "-----------------------------------------------------" << std::endl;
                read_nodes_coords_vector(nodes_coords_vector);
                std::cout << "-----------------------------------------------------" << std::endl;
                read_element_map(elem_nodes_vector);
                std::cout << "-----------------------------------------------------" << std::endl;
            }
            setup_distributed_mesh(nodes_coords_vector, elem_nodes_vector);
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
         * @return std::pair<NodeIdCoordsPairsVector, ElemIdNodeIdPairVector> the nodes_coords_vector and elem_nodes_vector of the line mesh. 
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

            NodeIdCoordsPairsVector nodes_coords_vector;
            ElemIdNodeIdPairVector elem_nodes_vector;
            nodes_coords_vector.reserve(divisions + 1);
            elem_nodes_vector.reserve(divisions);


            coords delta_xyz = (pts_coords[1] - pts_coords[0])/divisions;

            for (unsigned i = 0; i < divisions + 1; ++i)
            {
                nodes_coords_vector.push_back(std::make_pair(i + 1, pts_coords[0] + i*delta_xyz));
            }

            for (unsigned i = 0; i < divisions; ++i)
            {
                elem_nodes_vector.push_back(std::make_pair(i + 1, std::vector<unsigned>{i + 1, i + 2}));
            }     
            return std::make_pair(nodes_coords_vector, elem_nodes_vector);
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
         * @param sect a \ref BeamColumnFiberSection object that is used to initialise the beam-column elements.
         * 
        **/
       void create_distributed_line_mesh(int divisions, std::vector<coords> end_coords, ElementType elem_type, BeamColumnFiberSection& sect)
       {
            initialise_mpi_variables();
           element_type = elem_type;
           fiber_section = std::make_unique<BeamColumnFiberSection>(sect);
           std::pair<NodeIdCoordsPairsVector, ElemIdNodeIdPairVector> mesh_maps = map_a_line_mesh(divisions, end_coords);
           setup_distributed_mesh(mesh_maps.first, mesh_maps.second);
       }

        /**
         * 
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
         * @brief creates a line mesh map with a given number of divisions and end coordinates of the line. Does NOT take a gmsh mesh file. Mesh is distributed over a number of ranks.
         * @param elem_type an enum referring to the type of element that the mesh will include.
         * @param sect a \ref BasicSection object that is used to initialise the beam-column elements.
         * @param rank rank of the current process.
         * @param nranks number of ranks over which the mesh is distributed.
         * 
        **/
       void create_distributed_line_mesh(int divisions, std::vector<coords> end_coords, ElementType elem_type, BasicSection& sect)
       {
           initialise_mpi_variables();
           element_type = elem_type;
           basic_section = std::make_unique<BasicSection>(sect);
           std::pair<NodeIdCoordsPairsVector, ElemIdNodeIdPairVector> mesh_maps = map_a_line_mesh(divisions, end_coords);
           setup_distributed_mesh(mesh_maps.first, mesh_maps.second);
       }

        /**
         * @brief creates element objects following the \ref ElemIdNodeIdPairVector object format and adds them to \ref elem_vector.
         * 
         * @param nodes_coords_vector a vector of pairs mapping each node ids and coordinates.
         */
        void make_nodes (NodeIdCoordsPairsVector& nodes_coords_vector, bool interface = false)
        {
            if (!interface)
            {
                for (auto& node_data : nodes_coords_vector)
                {
                    node_vector.push_back(std::make_shared<Node>(node_data.first, node_data.second));
                }
            } else {
                for (auto& node_data : nodes_coords_vector)
                {
                    interface_node_vector.push_back(std::make_shared<Node>(node_data.first, node_data.second));
                }
            }
        }

        /**
         * @brief creates element objects following the \ref ElemIdNodeIdPairVector object format and adds them to \ref elem_vector.
         * 
         * @attention only one type of elements available now. Functionality limited to creating 2D beam-columns.
         * 
         * @param elem_nodes_vector a vector of pairs mapping elem ids and corresponding node ids.
         */
        void make_elements (ElemIdNodeIdPairVector& elem_nodes_vector)
        {
            #if VERBOSE
                std::cout << "Rank "<< rank << " -------------GlobalMesh::make_elements-------------" << std::endl;
            #endif
            std::vector<std::shared_ptr<Node>> elem_nodes;
            elem_nodes.reserve(2);
            for (auto& element_data : elem_nodes_vector)
            {
                elem_nodes.clear();
                for (auto& node_id: element_data.second)
                {
                    #if VERBOSE
                        std::cout << "Rank "<< rank << " - GlobalMesh::make_elements: searching for node id = " << node_id << " to make element " << element_data.first << std::endl;
                    #endif
                    std::shared_ptr<Node> node = get_node_by_record_id(node_id, "all");
                    elem_nodes.push_back(node);
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
         * @param nodes_coords_vector a vector of pairs mapping each node ids and coordinates.
         * @param elem_nodes_vector a vector of pairs mapping elem ids and corresponding 2 node ids.
         */
        void setup_mesh(NodeIdCoordsPairsVector nodes_coords_vector, ElemIdNodeIdPairVector elem_nodes_vector)
        {
            nnodes = nodes_coords_vector.size();
            nelems = elem_nodes_vector.size();
            node_vector.clear();
            node_vector.reserve(nnodes);
            elem_vector.clear();
            elem_vector.reserve(nelems);
            make_nodes(nodes_coords_vector);
            make_elements(elem_nodes_vector);
            sort_node_vector("all");
            count_dofs();            
        }

        void set_elem_type(ElementType elem_type)
        {
            element_type = elem_type;
        }
        void set_basic_section(BasicSection& sect)
        {
            basic_section = std::make_unique<BasicSection>(sect);
        }
        void set_fibre_section(BeamColumnFiberSection& sect)
        {
            fiber_section = std::make_unique<BeamColumnFiberSection>(sect);
        }


        /**
         * @brief populates a std::map object with a key being the node ID, and a value being the rank to which domain this node will belong using a naive algorithm where the nodes are divided amongst the ranks semi-equally with the last rank taking on any nodes that were not assigned already due to the remainder of nnodes/nranks. Also does the other way around populating a rank-node map.
         * 
         * @param node_rank_map std::map object with a key being the node ID, and a value being the rank to which domain this node will belong. Will be filled by reference.
         * @param node_id_set_owned_by_rank std::map object with a key being the rank and the value being a std::set of node IDs that belong to this rank. Will be filled by reference.
         * @param nodes_coords_vector the \ref nodes_coords_vector object which relates the node ID with their coordinates.
         * @param num_ranks the number of ranks over which the mesh is being decomposed.
         */
        void populate_node_rank_maps(std::map<unsigned, int>& node_rank_map, 
                                    std::set<unsigned>& node_id_set_owned_by_rank, 
                                    NodeIdCoordsPairsVector& nodes_coords_vector)
        {
            int nnodes = nodes_coords_vector.size();
            int nodes_per_rank = nnodes/num_ranks;

            for (int rank_i = 0; rank_i < num_ranks; ++rank_i)
            {
                if (rank_i != num_ranks - 1)
                {      
                    if (rank_i == rank)
                    {          
                        for (int i = rank_i*nodes_per_rank; i < (rank_i + 1)*nodes_per_rank; ++i)
                        {
                            node_rank_map[nodes_coords_vector[i].first] = rank_i;
                            node_id_set_owned_by_rank.insert(nodes_coords_vector[i].first);
                            
                        }
                    } else {
                        for (int i = rank_i*nodes_per_rank; i < (rank_i + 1)*nodes_per_rank; ++i)
                        {
                            node_rank_map[nodes_coords_vector[i].first] = rank_i;
                            
                        }
                    }
                }
                else
                {
                    if (rank_i == rank)
                    {
                    for (int i = rank_i*nodes_per_rank; i < nnodes; ++i)
                    {
                        node_rank_map[nodes_coords_vector[i].first] = rank_i;
                        node_id_set_owned_by_rank.insert(nodes_coords_vector[i].first);
                    }
                    } else {
                    for (int i = rank_i*nodes_per_rank; i < nnodes; ++i)
                    {
                        node_rank_map[nodes_coords_vector[i].first] = rank_i;
                    }
                    }
                }
            }
        }
        
        /**
         * @brief finds the IDs of all the nodes that need to be created on this rank. This includes nodes whose DoFs are owned by another rank but are needed for creating an interface element on the current rank. 
         * 
         * @param node_id_set_on_rank a set of unique node IDs that need to exist on the current rank.
         * @param elem_nodes_vector_on_rank a map whose key is element ID and value is the nodes of that element. Only elements that are created on this rank are included, including interface elements.
         */
        void find_rank_nodes(std::set<unsigned>& node_id_set_on_rank, 
                            ElemIdNodeIdPairVector& elem_nodes_vector_on_rank)
        {
            for (std::pair<unsigned,std::vector<unsigned>> elem_nodes_pair : elem_nodes_vector_on_rank)
            {
                for (unsigned node_id : elem_nodes_pair.second)
                {
                    node_id_set_on_rank.insert(node_id);
                }
            }  
        }

        /**
         * @brief 
         * 
         * @param interface_node_id_set_on_rank 
         * @param interface_elem_id_set_on_rank 
         * @param node_id_set_owned_by_rank 
         * @param elem_nodes_vector_on_rank 
         * @param rank 
         */
        void find_rank_interface_nodes_and_elems(std::set<unsigned>& interface_node_id_set_on_rank,
                                       std::set<unsigned>& interface_elem_id_set_on_rank,
                                       std::set<unsigned>& node_id_set_owned_by_rank,
                                       ElemIdNodeIdPairVector& elem_nodes_vector_on_rank)
        {
            for (std::pair<unsigned,std::vector<unsigned>> elem_nodes_pair : elem_nodes_vector_on_rank)
            {
                for (unsigned node_id : elem_nodes_pair.second)
                {
                    if (!node_id_set_owned_by_rank.count(node_id))
                    {
                        interface_node_id_set_on_rank.insert(node_id);
                        interface_elem_id_set_on_rank.insert(elem_nodes_pair.first);
                    }
                }
            }  
        }

        /**
         * @brief finds the node IDs that corresond to nodes that belong on this rank, but whose data will need to be communicated to other ranks. That is, nodes owned on current rank but are interface nodes on another.
         * @param neighbour_wanted_rank_node_id_map a std::map whose key is a neighbour's rank, and values are the node IDs that are owned by the current rank but will be an interface node on the neighbour.
         * @param interface_elem_id_set_on_rank a std::map whose key is a neighbour's rank, and values are the node IDs are owned by that neighbour and are interface nodes on this rank.
         * @param interface_node_id_set_on_rank 
         * @param elem_nodes_vector_on_rank 
         * @param node_rank_map 
         */
        void find_nodes_wanted_by_neighbours(std::map<int, std::set<unsigned>>& neighbour_wanted_rank_node_id_map,
                                             std::map<int, std::set<unsigned>>& wanted_from_neighbour_rank_node_id_map, 
                                             std::set<unsigned>& interface_elem_id_set_on_rank, 
                                             std::set<unsigned>& interface_node_id_set_on_rank, 
                                             ElemIdNodeIdPairVector& elem_nodes_vector_on_rank, 
                                             std::map<unsigned, int>& node_rank_map)
        {
            for (std::pair<unsigned,std::vector<unsigned>> elem_nodes_pair : elem_nodes_vector_on_rank)
            {
                if (interface_elem_id_set_on_rank.count(elem_nodes_pair.first))
                {
                    std::vector<int> neighbours_ranks;
                    for (unsigned node_id : elem_nodes_pair.second)
                    {
                        if (interface_node_id_set_on_rank.count(node_id))
                        {
                            neighbours_ranks.push_back(node_rank_map[node_id]);
                        } 
                    }
                    for (unsigned node_id : elem_nodes_pair.second)
                    {
                        int parent_rank = node_rank_map[node_id];
                        if (parent_rank == rank) // this is equivalent to the more intensive (!interface_node_id_set_on_rank.count(node_id)) as it is checking if the node is not an interface node which means it is owned by the current rank!
                        {    
                            for (int neighbor : neighbours_ranks)
                            {
                                #if (VERBOSE)
                                    std::cout << "rank "<< rank << " elem " << elem_nodes_pair.first << " node " << node_id << " with parent_rank " << parent_rank << " pushed to neighbour_wanted_rank_node_id_map[neighbor = " << neighbor<< "]." << std::endl;
                                #endif
                                neighbour_wanted_rank_node_id_map[neighbor].insert(node_id);
                            }
                        } else {
                            #if (VERBOSE)
                                std::cout << "rank " << rank << " elem " << elem_nodes_pair.first << " node " << node_id << " with parent_rank " << parent_rank << " pushed to wanted_from_neighbour_rank_node_id_map[parent_rank = " << parent_rank << "]." << std::endl;
                            #endif
                            wanted_from_neighbour_rank_node_id_map[parent_rank].insert(node_id);
                        }
                    }
                }
            } 
        }

        /**
         * @brief filters out \ref nodes_coords_vector to a \ref nodes_coords_vector_on_rank that corresponds only to the nodes on the current subdomain/rank. Does this by going over a set \ref node_id_set_on_rank that has a list of all the nodes on a given rank.
         * 
         * @param nodes_coords_vector_on_rank a vector of pairs containing node IDs and coordinates for each node on current rank.
         * @param nodes_coords_vector a vector of pairs containing node IDs and coordinates for each node in the global domain.
         * @param node_id_set_on_rank a std::set<unsigned> that has a nonrepeating set of node IDs for nodes that belong on current rank.
         * @param rank the rank that is calling this function.
         */
        void filter_node_vector(NodeIdCoordsPairsVector& nodes_coords_vector_on_rank,
                                NodeIdCoordsPairsVector& nodes_coords_vector, 
                                std::set<unsigned> node_id_set_on_rank)
        {
            for (unsigned node_id : node_id_set_on_rank)
            {
                auto node_coords_pair_it = std::find_if(std::begin(nodes_coords_vector), std::end(nodes_coords_vector), 
                                            [node_id](std::pair<unsigned, coords> node_coords_pair) 
                                            {
                                                return node_coords_pair.first == node_id;
                                            });
                nodes_coords_vector_on_rank.push_back(*node_coords_pair_it);
            }

        }
        /**
         * @brief checks if a given node is owned by this rank.
         * 
         * @param node_record_id the node record_ID that is being checked.
         * @return true if the calling rank owns this node.
         * @return false if the calling rank does not own this node.
         */
        bool owns_node_record_id(unsigned node_record_id)
        {
            return node_id_set_owned_by_rank.find(node_record_id) != node_id_set_owned_by_rank.end();
        }

        void read_nodal_U()
        {
            for (auto node: node_vector)
            {
                auto displacements = node->get_nodal_displacements();
                std::cout << "Rank " << rank << "-owned node " << node->get_id() << " has displacements: "; 
                print_container(displacements);
            }

            for (auto node: interface_node_vector)
            {
                auto displacements = node->get_nodal_displacements();
                std::cout << "Rank " << rank << " interface node " << node->get_id() << " has displacements: "; 
                print_container(displacements);
            }
        }

        /**
         * @brief populates the map node_element_map with connectivity information of the list of elements that are connected to each node. This is the inverse of \ref elem_nodes_vector which has pairs of element ID and the nodes that connect to the element.
         * 
         * @param node_element_map a map with node ID as key, and a set of elements that are connected to this node.
         * @param elem_nodes_vector a vector containing the nodal ids for each element in the domain.
         */
        void populate_node_element_map(std::map<unsigned, std::set<unsigned>>& node_element_map, 
                                        ElemIdNodeIdPairVector& elem_nodes_vector)
        {
            for (auto elem_nodes_pair : elem_nodes_vector)
            {
                for (auto node : elem_nodes_pair.second)
                {
                    node_element_map[node].insert(elem_nodes_pair.first);
                }
            }
        }

        /**
         * @brief finds the elements that belong on this current rank. Does this by checking the nodes that are owned by this rank, and for each node owned by this rank, finding which elements belong to it from the map node_element_set that was populated by \ref populate_node_element_map.
         * 
         * @param elem_id_set_on_rank a set of unique element IDs for elements that belong on this rank, including domain interface elements.
         * @param node_id_set_owned_by_rank a std::list which contents are the IDs of the nodes owned by this rank.
         * @param node_element_map a map which whose key is the node ID, and value is all the elements that connect to this node.
         */
        void find_rank_elements(std::set<unsigned>& elem_id_set_on_rank,
                                std::set<unsigned>&  node_id_set_owned_by_rank, 
                                std::map<unsigned, std::set<unsigned>> node_element_map)
        {
            for (unsigned node_id : node_id_set_owned_by_rank)
            {
                for (unsigned elem_id : node_element_map[node_id])
                {
                    elem_id_set_on_rank.insert(elem_id);
                }
            }   
        }

        /**
         * @brief filters out the pairs of element ID and element nodes for elements that reside on current rank from the global elem_nodes_vector.
         * 
         * @param elem_nodes_vector_on_rank vector of pairs of element ID and vector of connected node IDs for that will be created on the current rank.
         * @param elem_id_set_on_rank set of unique ID numbers for elements that need to be created on the current rank.
         * @param elem_nodes_vector vector of pairs of element ID and vector of connected node IDs for the global domain. 
         */
        void filter_element_vector(ElemIdNodeIdPairVector& elem_nodes_vector_on_rank,
                                std::set<unsigned>& elem_id_set_on_rank, 
                                ElemIdNodeIdPairVector& elem_nodes_vector)
        {
            for (unsigned elem_id : elem_id_set_on_rank)
            {

                auto elem_it = std::find_if(std::begin(elem_nodes_vector), std::end(elem_nodes_vector), 
                                        [elem_id](std::pair<unsigned, std::vector<unsigned>> elem_nodes_pair) 
                                        {
                                            return elem_nodes_pair.first == elem_id;
                                        });
                elem_nodes_vector_on_rank.push_back(*elem_it);
            }
        }

        /**
         * @brief finds the \ref max_num_stiffness_contributions by iterating over the nodes and finding the number of row contributions made by each of them.
         * 
         */
        void find_max_num_stiffness_contributions()
        {
            for (auto node_i : node_vector)
            {
                max_num_stiffness_contributions = std::max(max_num_stiffness_contributions, node_i->get_num_row_contributions());
            }
        }

        /**
         * @brief populates the global_mesh members: \ref nnodes, \ref nelems, \ref node_vector, and \ref elem_vector based on mesh (node and element) maps. Does the same for rank variants of these variables, and call \ref count_distributed_dofs to update the \ref nz_i of each node on the current rank.
         * 
         * @param nodes_coords_vector a vector of pairs mapping each node ids and coordinates.
         * @param elem_nodes_vector a vector of pairs mapping elem ids and corresponding 2 node ids.
         */
        void setup_distributed_mesh(NodeIdCoordsPairsVector& nodes_coords_vector, 
                                    ElemIdNodeIdPairVector& elem_nodes_vector)
        {
            nnodes = nodes_coords_vector.size();
            nelems = elem_nodes_vector.size();
            ranks_ndofs.resize(num_ranks);
            ranks_nnodes.resize(num_ranks);

            // sort nodes into the ranks that own them
            std::map<unsigned, int> node_rank_map; 
            node_id_set_owned_by_rank.clear();
            populate_node_rank_maps(node_rank_map, node_id_set_owned_by_rank, nodes_coords_vector);
            
            // create a map that links any nodes to the elements that connect to it.
            std::map<unsigned, std::set<unsigned>> node_element_map;
            populate_node_element_map(node_element_map, elem_nodes_vector);
            
            // Based on the nodes we currently have, find all elements that connect to any of these nodes.
            std::set<unsigned> elem_id_set_on_rank; 
            find_rank_elements(elem_id_set_on_rank, node_id_set_owned_by_rank, node_element_map);
            rank_nelems = elem_id_set_on_rank.size();
            // filter the elem_nodes_vector to only the members that belong on this rank (including those that are duplicated)
            ElemIdNodeIdPairVector elem_nodes_vector_on_rank;
            elem_nodes_vector_on_rank.reserve(rank_nelems);
            filter_element_vector(elem_nodes_vector_on_rank, elem_id_set_on_rank, elem_nodes_vector);

            // Based on the members that will be created on this rank, add the nodes that will also need to be created
            // interface and rank-owned node id set will be kept for BC and load operations handling.
            interface_node_id_set_on_rank.clear();
            std::set<unsigned> interface_elem_id_set_on_rank;
            find_rank_interface_nodes_and_elems(interface_node_id_set_on_rank,interface_elem_id_set_on_rank,node_id_set_owned_by_rank,elem_nodes_vector_on_rank);
            rank_nnodes = node_id_set_owned_by_rank.size();
            rank_interface_nnodes = interface_node_id_set_on_rank.size();

            // Find the node IDs that each rank may want from our nodes.

            find_nodes_wanted_by_neighbours(wanted_by_neighbour_rank_node_id_map, wanted_from_neighbour_rank_node_id_map, interface_elem_id_set_on_rank, interface_node_id_set_on_rank, elem_nodes_vector_on_rank, node_rank_map);
            
            // Add the nodes that officially belong to this rank to the rank nodes coords vector.
            NodeIdCoordsPairsVector nodes_coords_vector_on_rank;
            NodeIdCoordsPairsVector interface_nodes_coords_vector_on_rank;
            nodes_coords_vector_on_rank.reserve(rank_nnodes);
            interface_nodes_coords_vector_on_rank.reserve(rank_interface_nnodes);
            filter_node_vector(nodes_coords_vector_on_rank, nodes_coords_vector, node_id_set_owned_by_rank);
            filter_node_vector(interface_nodes_coords_vector_on_rank, nodes_coords_vector, interface_node_id_set_on_rank);

            // populate and sort the node and element object vectors for the rank
            node_vector.clear();
            node_vector.reserve(rank_nnodes);
            interface_node_vector.clear();
            interface_node_vector.reserve(rank_nnodes);

            elem_vector.clear();
            elem_vector.reserve(rank_nelems);
            make_nodes(nodes_coords_vector_on_rank);
            make_nodes(interface_nodes_coords_vector_on_rank, true);
            set_nodes_parent_ranks(rank, node_rank_map);

            if (VERBOSE)
            {
                std::cout << "Rank " << rank << " has created nodes: ";
                for (auto node : node_vector)
                    std::cout << node->get_id() << " ";
                std::cout << std::endl;
                std::cout << "Rank " << rank << " has created interface_nodes: ";
                for (auto node : interface_node_vector)
                    std::cout << node->get_id() << " ";
                std::cout << std::endl;
            }

            make_elements(elem_nodes_vector_on_rank);
            sort_node_vector("all");
            sort_element_vector();

            if (VERBOSE)
            {
                std::cout << "Rank " << rank << " has sorted its node vectors in anticipation of renumbering." << std::endl;
                std::cout << "Rank " << rank << " nodes are:";
                for (auto a_node : node_vector)
                {
                    std::cout << a_node->get_id() << " ";
                }
                std::cout << std::endl;
                std::cout << "Rank " << rank << " interface_nodes are:";
                for (auto a_node : interface_node_vector)
                {
                    std::cout << a_node->get_id() << " ";
                }
                std::cout << std::endl;
            }
            if (VERBOSE)
            {
                std::cout << "Rank" << rank << ": the functions renumber_nodes and sort_node_vector are next." << std::endl;
                std::cout << "--------------------------------------------------------------------------------------------------------" << std::endl;
                sleep(1);
            }     
            renumber_nodes();
            sort_node_vector("all");
            if (VERBOSE)
            {
                std::cout << "Rank " << rank << ": nodes renumbered and sorted. Exchanging interface node IDs next." << std::endl;
                std::cout << "--------------------------------------------------------------------------------------------------------" << std::endl;
                sleep(1);
                std::cout << "Rank " << rank << " nodes are:";
                for (auto a_node : node_vector)
                {
                    std::cout << a_node->get_id() << " ";
                }
                std::cout << std::endl;
                std::cout << "Rank " << rank << " interface_nodes are:";
                for (auto a_node : interface_node_vector)
                {
                    std::cout << a_node->get_id() << " ";
                }
                std::cout << std::endl;
            }
            if (VERBOSE)
            {
                std::cout << "Rank" << rank << ": the function exchange_interface_nodes_updated_ids is next." << std::endl;
                std::cout << "--------------------------------------------------------------------------------------------------------" << std::endl;
                sleep(1);
            }
            exchange_interface_nodes_updated_ids(wanted_by_neighbour_rank_node_id_map, wanted_from_neighbour_rank_node_id_map);
            // count the ndofs of each rank and assign each node an index that corresponds to the global matrices and vectors.
            if (VERBOSE)
            {
                std::cout << "Rank" << rank << ": the function count_distributed_dofs is next." << std::endl;
                std::cout << "--------------------------------------------------------------------------------------------------------" << std::endl;
                sleep(1);
            }
            // Is there a point in calling the count and exchange funciton here, or should it only be called by the `Model`?
            count_and_exchange_distributed_dofs();
        }

        /**
         * @brief counts the distributed DoFs and exchanges them updating the interface node DoFs as well.
         * 
         */
        void count_and_exchange_distributed_dofs()
        {
            count_distributed_dofs();
            if (VERBOSE)
            {
                std::cout << "Rank" << rank << ": the function exchange_interface_nodes_nz_i is next." << std::endl;
                std::cout << "--------------------------------------------------------------------------------------------------------" << std::endl;
                sleep(1);
            }
            exchange_interface_nodes_nz_i(wanted_by_neighbour_rank_node_id_map, wanted_from_neighbour_rank_node_id_map);
        }
        
        /**
         * @brief Filters a vector of node IDs to only include those matching the specified ownership type
         * 
         * @details Takes a vector of node IDs and returns a set containing only those IDs that exist in the 
         * specified node set (rank_owned nodes, interface nodes, or both). This docstring (not code) was AI generated but verified and modified.
         *
         * @param given_node_ids A vector of unsigned node IDs to be filtered
         * @param ownership A string specifying which nodes to match against:
         *                  "rank_owned" - Only nodes owned by this rank
         *                  "interface" - Only interface nodes on this rank
         *                  "all" - Both rank_owned and interface nodes
         * @return std::set<unsigned> A set containing the intersection of given_node_ids and the 
         *                            specified node set
         */
        template <typename Container>
        std::set<unsigned> filter_node_ids(Container& given_node_ids, std::string ownership = "rank_owned")
        {
            std::set<unsigned> given_id_set(given_node_ids.begin(), given_node_ids.end());
            std::set<unsigned> intersection;
            if (ownership == "rank_owned") {
                std::set_intersection(given_id_set.begin(), given_id_set.end(), 
                                    node_id_set_owned_by_rank.begin(), node_id_set_owned_by_rank.end(),
                                    std::inserter(intersection, intersection.begin()));
            }
            if (ownership == "interface") {
            std::set_intersection(given_id_set.begin(), given_id_set.end(), 
                                interface_node_id_set_on_rank.begin(), interface_node_id_set_on_rank.end(),
                                std::inserter(intersection, intersection.begin()));
            }
            if (ownership == "all")
            {
                std::set<unsigned> owned_intersection, interface_intersection;
                owned_intersection = filter_node_ids(given_node_ids, "rank_owned");
                interface_intersection = filter_node_ids(given_node_ids, "interface");
                intersection.insert(owned_intersection.begin(), owned_intersection.end());
                intersection.insert(interface_intersection.begin(), interface_intersection.end());
            }
            return intersection;
        }

        /**
         * @brief Set the nodes parent ranks based on current calling rank and the map node_rank_map.
         * 
         * @param rank rank that is currently calling this function.
         * @param node_rank_map std::map<unsigned, int> that links the node IDs with the ranks that own them.
         */
        void set_nodes_parent_ranks(int const rank, std::map<unsigned, int>& node_rank_map)
        {
            for (auto node: node_vector)
            {
                // nodes that are on the current rank have the same calling and parent ranks.
                node->set_parent_rank(rank, rank);
            }
            for (auto interface_node : interface_node_vector)
            {
                interface_node->set_parent_rank(node_rank_map[interface_node->get_record_id()], rank);
            }
        }
        /**
         * @brief exchnages the IDs of the interface nodes between neighbouring ranks.
         * 
         * @param wanted_by_neighbour_rank_node_id_map a std::map where the key is the neighbour rank and the value is a set of node IDs that this neighbour wants from the current rank
         * @param wanted_from_neighbour_rank_node_id_map  a std::map where the key is the neighbour rank and the value is a set of node IDs that this rank wants from the neighbour rank
         */
        void exchange_interface_nodes_updated_ids(std::map<int, std::set<unsigned>> wanted_by_neighbour_rank_node_id_map, 
                                                  std::map<int, std::set<unsigned>> wanted_from_neighbour_rank_node_id_map)
        {
            // get information about neighbours and buffer sizes
            int num_neighbours = 0;
            std::set<int> neighbours;
            std::map<int, int> neighbour_send_buffer_size_map;
            std::map<int, int> neighbour_receive_buffer_size_map;

            for (auto rank_nodes_set_pair : wanted_from_neighbour_rank_node_id_map)
            {
                neighbours.insert(rank_nodes_set_pair.first);
                neighbour_receive_buffer_size_map[rank_nodes_set_pair.first] = rank_nodes_set_pair.second.size();
            }

            for (auto rank_nodes_set_pair : wanted_by_neighbour_rank_node_id_map)
            {
                neighbours.insert(rank_nodes_set_pair.first);
                neighbour_send_buffer_size_map[rank_nodes_set_pair.first] = rank_nodes_set_pair.second.size();
            }
            num_neighbours = neighbours.size();
            if (VERBOSE)
            {
                std::cout << "rank " << rank << " has " << num_neighbours << " neighbours which are: ";
                print_container(neighbours);
                for (auto neighbour_i : neighbours)
                {
                    std::cout << "Rank " << rank << " - wanted_from_neighbour_rank_node_id_map[" << neighbour_i << "] = ";
                    print_container(wanted_from_neighbour_rank_node_id_map[neighbour_i]);

                    std::cout << "Rank " << rank << " - wanted_by_neighbour_rank_node_id_map[" << neighbour_i << "] = ";
                    print_container(wanted_by_neighbour_rank_node_id_map[neighbour_i]);   
                }
                sleep(1);
            }

            // setup send buffers
            std::map<int, std::vector<unsigned>> rank_ids_send_buffers_map;
            for (auto rank_buffer_size_pair : neighbour_send_buffer_size_map)
            {
                rank_ids_send_buffers_map[rank_buffer_size_pair.first].reserve(rank_buffer_size_pair.second);
            }

            // setup receive buffers
            std::map<int, std::vector<unsigned>> rank_id_receive_buffers_map;
            for (auto rank_buffer_size_pair : neighbour_receive_buffer_size_map)
            {
                rank_id_receive_buffers_map[rank_buffer_size_pair.first].resize(rank_buffer_size_pair.second);
            }
            

            
            if (VERBOSE)
                std::cout << "rank " << rank << " populating send buffers started." << std::endl;
            // populate send buffers --- std::map<int, std::set<unsigned>> wanted_by_neighbour_rank_node_id_map
            for (auto rank_node_set : wanted_by_neighbour_rank_node_id_map)
            {
                if (VERBOSE)
                    std::cout << "rank " << rank << " looping over set for neighbor " << rank_node_set.first << std::endl;
                for (auto node_id : rank_node_set.second)
                {
                    if (VERBOSE)
                        std::cout << "rank " << rank << " looking for node with record_id " << node_id << " for neighbor " << rank_node_set.first << std::endl;
                    unsigned renumbered_id = get_node_by_record_id(node_id, "rank_owned")->get_id();
                    rank_ids_send_buffers_map[rank_node_set.first].push_back(renumbered_id);
                }
            }
            if (VERBOSE)
                std::cout << "rank " << rank << " populating send buffers completed." << std::endl;

            if (VERBOSE)
                std::cout << "rank " << rank << " exchanging nodes ids with MPI started." << std::endl;
            // exchange node ids
            #ifdef WITH_MPI
             for (int neighbor_rank : neighbours)
            {

                auto receive_buffer = rank_id_receive_buffers_map[neighbor_rank].data();
                auto send_buffer = rank_ids_send_buffers_map[neighbor_rank].data();
                int send_buffer_size = neighbour_send_buffer_size_map[neighbor_rank];
                int receive_buffer_size = neighbour_receive_buffer_size_map[neighbor_rank];
                if (VERBOSE)
                {
                    std::cout << "GlobalMesh::exchange_interface_nodes_updated_ids: Rank " << rank << " (sending, receiving) = (" << send_buffer_size << ", " << receive_buffer_size << ") objects to neighbor_rank " << neighbor_rank <<
                        "; send buffer size is " << rank_ids_send_buffers_map[neighbor_rank].size() <<
                        " and receive_buffer size is " << rank_id_receive_buffers_map[neighbor_rank].size() << "." << std::endl;
                }
                MPI_Sendrecv(send_buffer, send_buffer_size, MPI_UNSIGNED, neighbor_rank, 0,
                            receive_buffer, receive_buffer_size, MPI_UNSIGNED, neighbor_rank, 0, 
                            MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            
            #endif

            
            if (VERBOSE)
                std::cout << "rank " << rank << " exchanging nodes ids with MPI completed." << std::endl;
            if (VERBOSE)
                std::cout << "rank " << rank << " copying back from receive buffer started." << std::endl;
            // copy back from the receive buffer into the nodes new ids
            for (auto rank_id_buffer_pair : rank_id_receive_buffers_map)
            {
                auto node_id_set_iterator = wanted_from_neighbour_rank_node_id_map[rank_id_buffer_pair.first].begin();
                for (int i = 0; i < rank_id_buffer_pair.second.size(); ++i)
                {
                    unsigned new_id = rank_id_buffer_pair.second[i];
                    unsigned original_record_id = *node_id_set_iterator;
                    if (VERBOSE)
                     std::cout << "rank " << rank << " looking for node " << original_record_id << " to assign it new id = " << new_id << std::endl;
                    
                     get_node_by_record_id(original_record_id, "interface")->set_id(new_id);
                    ++node_id_set_iterator;
                }
            }
            if (VERBOSE)
            {
                std::cout << "rank " << rank << " exchange_interface_nodes_updated_ids all completed============================." << std::endl;
                sleep(1);
            }
                
        }


        void read_node_ids(bool interface = false)
        {
            if (!interface)
            {
            for (auto a_node : node_vector)
            {
                std::cout << "read_node_ids: Rank " << rank << " has node ID = " << a_node->get_id() << " and record_id = " << a_node->get_record_id() << std::endl;
            }
            } else {
            for (auto a_node : interface_node_vector)
            {
                std::cout << "read_node_ids: Rank " << rank << " has interface node ID = " << a_node->get_id() << " and record_id = " << a_node->get_record_id() << std::endl;
            }
        }
        }

        std::vector<int> get_ranks_ndofs_vector() {return ranks_ndofs;}

        // void print_local_nodes()
        // {
        //     std::cout << 
        // }
        // void print_local_ndofs()
        // {

        // }


        /**
         * @brief exchnages the nz_i of the interface nodes between neighbouring ranks.
         * 
         * @param wanted_by_neighbour_rank_node_id_map 
         * @param wanted_from_neighbour_rank_node_id_map 
         */
        void exchange_interface_nodes_nz_i(std::map<int, std::set<unsigned>> wanted_by_neighbour_rank_node_id_map, 
                                                  std::map<int, std::set<unsigned>> wanted_from_neighbour_rank_node_id_map)
        {
            // get information about neighbours and buffer sizes
            int num_neighbours = 0;
            std::set<int> neighbours;
            std::map<int, int> neighbour_send_buffer_size_map;
            std::map<int, int> neighbour_receive_buffer_size_map;

            for (auto rank_nodes_set_pair : wanted_from_neighbour_rank_node_id_map)
            {
                neighbours.insert(rank_nodes_set_pair.first);
                neighbour_receive_buffer_size_map[rank_nodes_set_pair.first] = rank_nodes_set_pair.second.size();
            }
            for (auto rank_nodes_set_pair : wanted_by_neighbour_rank_node_id_map)
            {
                neighbours.insert(rank_nodes_set_pair.first);
                neighbour_send_buffer_size_map[rank_nodes_set_pair.first] = rank_nodes_set_pair.second.size();
            }

            num_neighbours = neighbours.size();


            // setup send buffers
            std::map<int, std::vector<int>> rank_nz_i_send_buffers_map;
            for (auto rank_buffer_size_pair : neighbour_send_buffer_size_map)
            {
                rank_nz_i_send_buffers_map[rank_buffer_size_pair.first].reserve(rank_buffer_size_pair.second);
            }

            // setup receive buffers
            std::map<int, std::vector<int>> rank_nz_i_receive_buffers_map;
            for (auto rank_buffer_size_pair : neighbour_receive_buffer_size_map)
            {
                rank_nz_i_receive_buffers_map[rank_buffer_size_pair.first].resize(rank_buffer_size_pair.second);
            }


            // populate send buffers --- std::map<int, std::set<unsigned>> wanted_by_neighbour_rank_node_id_map
            for (auto rank_node_set : wanted_by_neighbour_rank_node_id_map)
            {
                for (auto node_id : rank_node_set.second)
                {
                    unsigned node_nz_i = get_node_by_record_id(node_id, "rank_owned")->get_nz_i();
                    rank_nz_i_send_buffers_map[rank_node_set.first].push_back(node_nz_i);
                }
            }

            // exchange node ids
            #ifdef WITH_MPI
             for (int neighbor_rank : neighbours)
            {
                auto receive_buffer = rank_nz_i_receive_buffers_map[neighbor_rank].data();
                auto send_buffer = rank_nz_i_send_buffers_map[neighbor_rank].data();
                int send_buffer_size = neighbour_send_buffer_size_map[neighbor_rank];
                int receive_buffer_size = neighbour_receive_buffer_size_map[neighbor_rank];
                if (VERBOSE)
                {
                    std::cout << "GlobalMesh::exchange_interface_nodes_nz_i: Rank " << rank << " (sending, receiving) = (" << send_buffer_size << ", " << receive_buffer_size << ") objects to neighbor_rank " 
                        << neighbor_rank << "; send buffer size is " << rank_nz_i_send_buffers_map[neighbor_rank].size() <<
                        " and receive_buffer size is " << rank_nz_i_receive_buffers_map[neighbor_rank].size() << "." << std::endl;
                }
                MPI_Sendrecv(send_buffer, send_buffer_size, MPI_INT, neighbor_rank, 0,
                            receive_buffer, receive_buffer_size, MPI_INT, neighbor_rank, 0, 
                            MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            #endif
            // copy back from the receive buffer into the nodes new ids
            for (auto rank_id_buffer_pair : rank_nz_i_receive_buffers_map)
            {
                auto node_id_set_iterator = wanted_from_neighbour_rank_node_id_map[rank_id_buffer_pair.first].begin();
                for (int i = 0; i < rank_id_buffer_pair.second.size(); ++i)
                {
                    int new_nz_i = rank_id_buffer_pair.second[i];
                    unsigned original_record_id = *node_id_set_iterator;
                    get_node_by_record_id(original_record_id, "interface")->set_nz_i(new_nz_i);
                    ++node_id_set_iterator;
                }
            }
        }

        /**
         * @brief renumbers the nodes on the current rank so that they are continuous and take into account the number of nodes on all other ranks as well. Requires an MPI_Allgather operation, and requires communication to update the node number for nodes that do not belong on the current rank.
         * 
         */
        void renumber_nodes()
        {
            unsigned* ranks_nnodes_ptr = ranks_nnodes.data();
            #ifdef WITH_MPI
            if (VERBOSE)
            {
                std::cout << "GlobalMesh::renumber_nodes: Rank " << rank << " sending via MPI_Allgather " << 1 << " objects to all ranks;" 
                    << " send buffer size is " << 1 <<
                    " and receive_buffer size is: " << ranks_nnodes.size() << "." << std::endl;
            }
            // Find out what the rank_nnodes is for each rank
            // MPI_Barrier(MPI_COMM_WORLD);
            MPI_Allgather(&rank_nnodes, 1, MPI_UNSIGNED,
                        ranks_nnodes_ptr, 1, MPI_UNSIGNED, 
                        MPI_COMM_WORLD);
            #endif

            // update the node id for each node
            // Step 1: find the number to udpate the rank_starting_node_id
            rank_starting_node_id = 1;
            if (rank > 0)
            {
                for (int i = 0; i < rank; ++i)
                {
                    rank_starting_node_id += ranks_nnodes[i];
                }
            }
            if (VERBOSE)
            {
                std::cout << "GlobalMesh::renumber_nodes: Rank " << rank << ": rank_starting_node_id is: "<< rank_starting_node_id << "." << std::endl;
            }
            // Step 2: loop over the nodes and update them.
            unsigned temp_id = rank_starting_node_id;
            for (auto& node: node_vector)
            {
                if (VERBOSE)
                {
                    std::cout << "Rank " << rank << ": Node " << node->get_id() << " will be renumbered to " << temp_id << std::endl;
                }
                node->set_id(temp_id);
                ++temp_id;
                if (VERBOSE)
                {
                    std::cout << "Rank " << rank << ": Node " << node->get_id() << " has record_id = " << node->get_record_id() << std::endl;
                }
            }

        }
        /**
         * @brief counts the active DOFs in the mesh by going over all the nodes and getting the number of active freedoms. This is done over a distributed domain by first carrying out the operation locally, then getting the number of DoFs on each rank via an MPI_Allgather call. The \ref rank_starting_nz_i is udpated for the current rank by summing all ndofs of the ranks lower than it. This rank_starting_nz_i is then used to update the nz_i of each node which was initialised with a local nz_i.
         */
        void count_distributed_dofs()
        {
            int* ranks_ndofs_ptr = ranks_ndofs.data();
            rank_ndofs = 0;
            for (auto& node: node_vector)
            {
                node->set_nz_i(rank_ndofs);
                if (VERBOSE)
                {
                    std::cout << "Node " << node->get_id() << " has nz_i = " << node->get_nz_i() << std::endl;
                }
                rank_ndofs += node->get_ndof();
            }
            #ifdef WITH_MPI
            // Find out what the rank_ndofs is for each rank
            if (VERBOSE)
            {
                std::cout << "GlobalMesh::count_distributed_dofs: Rank " << rank << " sending via MPI_Allgather " << 1 << " objects to all ranks;" 
                    << " send buffer size is " << 1 <<
                    " and receive_buffer size is: " << ranks_ndofs.size() << "." << std::endl;
            }
            MPI_Allgather(&rank_ndofs, 1, MPI_INT,
                        ranks_ndofs_ptr, 1, MPI_INT, 
                        MPI_COMM_WORLD);
            #endif

            // update the nz_it for each node
            // Step 1: find the number to udpate the nz_i of the nodes
            rank_starting_nz_i = 0;
            if (rank > 0)
            {
                for (int i = 0; i < rank; ++i)
                {
                    rank_starting_nz_i += ranks_ndofs[i];
                }
            }
            // Step 2: loop over the nodes and update them.
            for (auto& node: node_vector)
            {
                node->increment_nz_i(rank_starting_nz_i);
                if (VERBOSE)
                {
                    std::cout << "Node " << node->get_id() << " has nz_i = " << node->get_nz_i() << std::endl;
                }
            }

            // Step 3: count the total ndofs of the whole problem
            ndofs = 0;
            for (int rank_i = 0; rank_i < num_ranks; ++rank_i)
            {
                ndofs += ranks_ndofs[rank_i];
            }
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
                    std::cout << "Node " << node->get_id() << " has nz_i = " << node->get_nz_i() << std::endl;
                }
                ndofs += node->get_ndof();
            }
        }

        int get_ndofs() {return ndofs;}
        int get_rank_ndofs() {return rank_ndofs;}
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
        void fix_node(int const id, int const dof)
        {
            auto node_ptr = get_node_by_record_id(id, "all");
            if (dof < 0)
            {
                std::cout << "Fixing all DoFs of node " << id << std::endl;
                node_ptr->fix_all_dofs();
                node_ptr->print_inactive_dofs();
            } else {
                node_ptr->fix_dof(dof);
            }    
        }
        /**
         * @brief calls a node's \ref Node::add_nodal_load command and loads it.
         * 
         */
        void load_node(int id, int dof, real load)
        {
            get_node_by_record_id(id, "rank_owned")->add_nodal_load(load, dof);
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
            get_node_by_record_id(id, "rank_owned")->increment_nodal_load(dP, dof);
        }

        /**
         * @brief keeps track of the history of a nodal DoF.
         * 
         * @param id 
         * @param dof 
         * @param history vector where nodal DoF history is tracked.
         */
        void track_nodal_dof(int const id, int const dof, std::vector<real>& history)
        {
            auto node_ptr = get_node_by_record_id(id, "rank_owned");
            std::array<real, 6> nodal_displacements = node_ptr->get_nodal_displacements();
            history.push_back(nodal_displacements[dof]);
        }

        /**
         * @brief get node shared_ptr by id.
         * @param id id of the node to get from the \ref node_vector or \ref interface_node_vector or both.
         * @param search_target a std::string that is either "all", "rank_owned", or "interface" that specifies where to search for the nodes.
         * @return std::shared_ptr<Node> a shared pointer to the node with the given id.
         */
        std::shared_ptr<Node> get_node_by_id(int id, std::string search_target = "all")
        {
            std::vector<std::shared_ptr<Node>>::iterator node_it;
            if (search_target == "rank_owned")
            {
                 node_it = get_record_id_iterator<std::vector<std::shared_ptr<Node>>::iterator, std::vector<std::shared_ptr<Node>>>(id, node_vector);
            }
            if (search_target == "interface")
            {
                 node_it = get_record_id_iterator<std::vector<std::shared_ptr<Node>>::iterator, std::vector<std::shared_ptr<Node>>>(id, interface_node_vector);
            }
 
            if (search_target == "all")
            {
                 node_it = get_record_id_iterator<std::vector<std::shared_ptr<Node>>::iterator, std::vector<std::shared_ptr<Node>>>(id, node_vector);
 
                 if (node_it == node_vector.end())
                 {
                     node_it = get_record_id_iterator<std::vector<std::shared_ptr<Node>>::iterator, std::vector<std::shared_ptr<Node>>>(id, interface_node_vector);
                 }
            }
            if ((node_it != node_vector.end()) && (node_it != interface_node_vector.end()))
            {
             return *node_it;
            } else {
             std::cout << "Could not find node with record_id " << id << " in " << search_target << " node vectors." << std::endl;
             exit(1);
            }
        }

        /**
         * @brief get node shared_ptr by record_id.
         * @param record_id id of the node to get from the \ref node_vector or \ref interface_node_vector or both.
         * @param search_target a std::string that is either "all", "rank_owned", or "interface" that specifies where to search for the nodes.
         * @return std::shared_ptr<Node> a shared pointer to the node with the given record_id.
         */
        std::shared_ptr<Node> get_node_by_record_id(int record_id, std::string search_target = "all")
        {
            std::vector<std::shared_ptr<Node>>::iterator node_it;
           if (search_target == "rank_owned")
           {
                node_it = get_record_id_iterator<std::vector<std::shared_ptr<Node>>::iterator, std::vector<std::shared_ptr<Node>>>(record_id, node_vector);
                if (node_it != node_vector.end())
                {
                    return *node_it;
                }
            }
           if (search_target == "interface")
           {
                if (VERBOSE)
                {
                    std::cout << "get_node_by_record_id::Rank " << rank << ": is looking for record_id " << record_id << " in " << search_target << std::endl; 
                    read_node_ids(true);
                }
                
                node_it = get_record_id_iterator<std::vector<std::shared_ptr<Node>>::iterator, std::vector<std::shared_ptr<Node>>>(record_id, interface_node_vector);
                if (node_it != interface_node_vector.end())
                {
                    return *node_it;
                }
            }

           if (search_target == "all")
           {
                node_it = get_record_id_iterator<std::vector<std::shared_ptr<Node>>::iterator, std::vector<std::shared_ptr<Node>>>(record_id, node_vector);

                if (node_it != node_vector.end())
                {
                    return *node_it;
                } else {
                    node_it = get_record_id_iterator<std::vector<std::shared_ptr<Node>>::iterator, std::vector<std::shared_ptr<Node>>>(record_id, interface_node_vector);
                    if (node_it != interface_node_vector.end())
                    {
                        return *node_it;
                    }
                }
           }
     
            std::cout << "Rank " << rank << ": Could not find node with record_id " << record_id << " in " << search_target << " node vectors." << std::endl;
            exit(1);
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
         * @brief makes each node calculate how much it will contribute to \f$\boldsymbol{P}\f$ by calling \ref Node::compute_global_load_triplets().
         * 
         */
        void calc_nodal_contributions_to_P()
        {
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
            #ifdef KOKKOS
                  Kokkos::parallel_for( "GlobalMesh::update_elements_states", elem_vector.size(), KOKKOS_LAMBDA (int i) {
                        elem_vector[i]->update_state();
                        elem_vector[i]->calc_global_stiffness_triplets();
                    });
            #else
                #pragma omp parallel for
                for (auto& elem: elem_vector)
                {
                    elem->update_state();
                    elem->calc_global_stiffness_triplets();
                }
            #endif

        }

        void update_element_sections_starting_states()
        {
            #ifdef KOKKOS
                Kokkos::parallel_for( "GlobalMesh::update_element_sections_starting_states", elem_vector.size(), KOKKOS_LAMBDA (int i) {
                    elem_vector[i]->update_section_starting_state();
                });
            #else
                #pragma omp parallel for
                for (auto& elem: elem_vector)
                {
                    elem->update_section_starting_state();
                }
            #endif
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
        

        /**
         * @brief sorts \ref node_vector and/or \ref interface_node_vector by ID.
         * @param vector_to_sort a std::string that is either "all", "rank_owned", or "interface" that specifies which vector of nodes to sort.
         * @warning depends on defining the "<" operator in \ref Node.
         */
        void sort_node_vector(std::string vector_to_sort="all")
        {
            if (vector_to_sort == "all")
            {
                sort_node_vector("rank_owned");
                sort_node_vector("interface");
            }
            else if (vector_to_sort == "rank_owned")
            {
                std::sort(node_vector.begin(), node_vector.end(), 
                        [](const auto& a, const auto& b) 
                        {
                        return *a < *b;
                        }
                        );
            }
            else if (vector_to_sort == "interface")
            {
                std::sort(interface_node_vector.begin(), interface_node_vector.end(), 
                        [](const auto& a, const auto& b) 
                        {
                        return *a < *b;
                        }
                        );
            }
            else 
            {
                std::cout << "ERROR: GlobalMesh::sort_node_vector can only accept all, rank_owned, or interface as vector_to_sort. Received " << vector_to_sort << std::endl; 
                exit(1);
            }
        }
        /**
         * @brief sorts \ref elem_vector by ID.
         */
        void sort_element_vector()
        {
        std::sort(elem_vector.begin(), elem_vector.end(), 
                [](const auto& a, const auto& b) 
                {
                return a->get_id() < b->get_id();
                }
                );
        }
        int const get_num_elems() const {return nelems;}
        int const get_num_nodes() const {return nnodes;}
        int const get_rank_num_nodes() const {return rank_nnodes;}
        int const get_rank_num_elems() const {return rank_nelems;}

        int const count_nodes_vector() const {return node_vector.size();}
        int const count_interface_nodes_vector() const {return interface_node_vector.size();}
        int const count_elem_vector() const {return elem_vector.size();}
        /**
         * @brief checks if \ref node_vector and/or \ref interface_node_vector contain a given set of node ids.alignas
         * 
         * @warning if there are fewer node_ids than in the mesh vectors and all of them are within the mesh, then this returns true although there are IDs in the mesh that are not in the passed node_ids set. 
         * 
         * @param node_ids a std::set of node IDs to check whether or not they are created in the mesh.          
         * @param search_target a std::string that is either "all", "rank_owned", or "interface" that specifies where to search for the nodes.
         * @return true if the node_ids are within the node vectors.
         * @return false if one or more IDs from node_ids are not within the node vectors.
         */
        bool contains_nodes(std::set<unsigned> node_ids, std::string search_target) 
        {
            bool contains_nodes;
            if (search_target == "all" || search_target == "rank_owned")
            {
                contains_nodes = contains_id(node_ids, node_vector);
                if (!contains_nodes && search_target == "all")
                {
                    contains_nodes = contains_id(node_ids, interface_node_vector);
                }
            }
            if (search_target == "interface")
            {
                contains_nodes = contains_id(node_ids, interface_node_vector);
            }
            return contains_nodes;
        }

        bool contains_elements(std::set<unsigned> elem_ids) 
        {
            return contains_id(elem_ids, elem_vector);
        }

        int get_rank_ndofs() const {return rank_ndofs;}
        int get_rank_starting_nz_i() const {return rank_starting_nz_i;}
};
#endif
