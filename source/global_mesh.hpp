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

/**
 * @brief std vector of pairs each of which has an id and a 3-item coords vector
 * 
 */
using gmsh_node_map = std::vector<std::pair<size_t, coords>>;

/**
 * @brief std vector of pairs each of which has an id and a vector of node ids related to the element id
 * 
 */
using gmsh_elem_map = std::vector<std::pair<size_t, std::vector<size_t>>>;


/**
 * @brief contains the objects and functions operating on the global stiffness, and displacement and force vectors. Also contains the references to the nodes and elements.
 * 
 * @attention does not use the global matrices/vectors, so they should be removed as members
 */
class GlobalMesh {
    private: 
        int nnodes = 0; /**< number of nodes in the mesh*/
        int ndofs = 0; /**< number of DOFs in the mesh*/
        int nelems = 0; /**< number of elements in the mesh*/
        std::vector<std::shared_ptr<Node>> node_vector;  /**< a vector of shared ptrs referring to all the nodes in the problem*/
        std::vector<std::shared_ptr<Basic2DBeamElement>> elem_vector; /**< a vector of shared ptrs referring to all the elements in the problem*/

        spmat K;
        spvec P;
        vec U; 

    public:
        friend class Assembler;
        void open_mesh_file(std::string const mesh_file);
        /**
         * @brief uses gmsh API to read the nodes and populate a map with ids and coordinates
         * 
         * @return gmsh_node_map a vector of pairs mapping each node ids and coordinates
         */
        gmsh_node_map read_nodes();

        /**
         * @brief reads the elements from the mesh file and populates the map linking element id with its node ids
         * 
         * @return gmsh_elem_map elem_map a vectpr of pairs mapping elem ids and corresponding node ids
         */
        gmsh_elem_map read_elements();

        /**
         * @brief creates element objects following the \ref gmsh_elem_map object format and adds them to \ref elem_vector
         * 
         * @param node_map a vector of pairs mapping each node ids and coordinates
         */
        void make_nodes (gmsh_node_map node_map);

        /**
         * @brief creates element objects following the \ref gmsh_elem_map object format and adds them to \ref elem_vector
         * 
         * @attention only one type of elements available now. Functionality limited to creating 2D beam-columns
         * 
         * @param elem_map a vector of pairs mapping elem ids and corresponding node ids
         */
        void make_elements (gmsh_elem_map elem_map);
        void close_mesh_file();

        /**
         * @brief populates the global_mesh members: \ref nnodes, \ref nelems, \ref node_vector, and \ref elem_vector
         * 
         * @param mesh_file a string that is the mesh file name
         */
        void setup_mesh(std::string const mesh_file);

        /**
         * @brief counts the active DOFs in the mesh by going over all the nodes and getting the number of active freedoms
         * 
         */
        void count_dofs();
        void print_info();

        void fix_node(int id, int dof);
        void calc_global_contributions() {
            std::cout << "Calc_global_contirbutions: There are " << ndofs << " active DoFs in the mesh." << std::endl;
            for (auto elem: elem_vector) 
            {   
                elem->map_stiffness();
                elem->calc_K_global();
            }
        }
        void solve_for_U();
        int const get_num_elems() const {return nelems;}
};

/**
 * @brief Get iterator for a node or element by searching for their id
 * 
 * @details searches by checking if the id of the element corrsponds closely to its location
 * in the vector. if not, it goes either up or down to keep checking.
 * this "search" is inefficient compared to other search algorithms such as std::find_if
 * or std::lower_bound for the general case, but is more efficient considering the average case we 
 * actually care about: a sorted vector of nodes/elems that is almost always continguous
 * 
 * @attention there might be a potential infinite loop in this search
 * 
 * @todo create test cases that searches when there is no node or element with the given id
 * in the searched container
 * 
 * @tparam Iterator stl-compatible iterator corresponding to the stl compatible container
 * @tparam Container any container with stl compatible interface
 * @param id unique id of node or vector to search
 * @param a_vec the container containing the nodes or the elements
 * @return Iterator 
 */
template <typename Iterator, typename Container>
Iterator get_id_iterator(int id, Container& a_vec)
{
    auto itr = std::begin(a_vec) + (id - 1);
    int check_id = ((*itr)->get_id());
    if (check_id > id)
    {
        while (check_id != id && (itr > std::begin(a_vec)))
        {
            --itr;
            check_id = ((*itr) -> get_id());
        }
    } else if (check_id < id) {
        while (check_id != id && (itr < std::end(a_vec)))
        {
            ++itr;
            check_id = ((*itr) -> get_id());
        }
    }
    if (check_id == id)
    {
        return itr;
    } else 
    {
        std::cout << "could not find item with id " << id << " in vector." << std::endl;
        std::exit(1);
    }
    
}
/**
 * @brief helps assemble the global matrices
 * 
 * @details supposed to help with organising the code by separating assembly; right now it is not helpful and only contains one big function and redundunt member variables that appeared in \ref GlobalMesh
 * 
 */
class Assembler {
    private:
        spmat K;
        spvec P;
        // TODO: Figure out if U should be sparse or dense!
        // vec U;
        spvec U;
    public:
        friend class BasicSolver;
        /**
         * @brief retrives global contributions from all elements
         * 
         * @details creates triplets by retrieving all global contributions from
         * the elements and then uses the triplets to creat the global sparse stiffness
         * matrix. Also allocates the force and displacement vectors.
         * 
         * @attention sets a force of -1e4 in one of the locations.
         * 
         * @todo add a function to read and apply forces to nodes
         * 
         * @param glob_mesh takes the global_mesh object as input to get the counters and containers for nodes and elements.
         */
        void assemble_global_contributions(GlobalMesh& glob_mesh);


};

/**
 * @brief place-holder class for solvers
 * 
 */
class BasicSolver {
    public:
        /**
         * @brief solves for U using the global matrices contained in \ref Assembler; uses Eigen's SparseLU solver
         * 
         * @param assembler 
         */
        void solve_for_U(Assembler& assembler);
};
bool check_matrix(spmat A);
bool has_zero_row(spmat A);