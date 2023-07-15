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

using gmsh_node_map = std::vector<std::pair<size_t, coords>>;
using gmsh_elem_map = std::vector<std::pair<size_t, std::vector<size_t>>>;

class global_mesh {
    private: 
        int nnodes = 0;
        int ndofs = 0;
        int nelems = 0;
        std::vector<std::shared_ptr<Node>> node_vector;
        std::vector<std::shared_ptr<Basic2DBeamElement>> elem_vector;

        spmat K;
        spvec P;
        vec U; 

    public:
        void open_mesh_file(std::string const mesh_file);
        gmsh_node_map read_nodes();
        gmsh_elem_map read_elements();
        void make_nodes (gmsh_node_map node_map);
        void make_elements (gmsh_elem_map elem_map);
        void close_mesh_file();
        void setup_mesh(std::string const mesh_file);
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
        void assemble_global_contributions();
        void solve_for_U();
        int const get_num_elems() const {return nelems;}
};

template <typename Iterator, typename Container>
Iterator get_id_iterator(int id, Container& a_vec)
{
    auto itr = std::begin(a_vec) + (id - 1);
    int check_id = ((*itr)->get_id());
    // this "search" is inefficient compared to other search
    // algorithms such as std::find_if or std::lower_bound
    // for the general case, but is more efficient considering
    // the average case we actually care about: a sorted
    // vector of nodes that is almost always continguous
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
bool check_matrix(spmat A);
bool has_zero_row(spmat A);

