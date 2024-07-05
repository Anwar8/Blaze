/**
 * @file BeamElementCommonInterface.hpp
 * @brief includes the common interfaces that are to be used by ALL derived beam-column element classes.
 */

#ifndef BEAM_ELEMENT_COMMON_INTERFACE_HPP
#define BEAM_ELEMENT_COMMON_INTERFACE_HPP

#include "BeamElementBaseClass.hpp"

/**
 * @brief Common interface for all beam-column elements.
 * 
 * @details This class has implementations for parts of \ref BeamElementBaseClass interface that are common across all beam-column elements. It is a virtual class that is to be inherited by all beam-column elements, and is not meant to be instantiated on its own. 
 * That being said, it does not contain the common variables that are used by all beam-column elements, only implementations for some functions (particularly related to mapping across local and global DOFs). All variables are in \ref BeamElementBaseClass.
 * 
 */
class BeamElementCommonInterface : public BeamElementBaseClass {
    protected:
        /**
         * @name basic_information
         * @brief the basic data about the generic beam element.
         * all basic_information is inherited from \ref BeamElementBaseClass.
         */
        //@{      
        //@}

        /**
         * @name beam_basic_objects
         * @brief basic objects needed by the beam-column elements. Section, shape function, transformation, etc.
         * all beam_basic_objects are inherited from \ref BeamElementBaseClass
         */
        //@{
        //@}

        /**
         * @name beam_state_containers
         * @brief the containers for the beam state such as displacement, strain, force, etc.
         * all beam_state_containers are inherited from \ref BeamElementBaseClass.
         */
        //@{
        //@}

    public:
    /**
     * @name element initialisers
     * @brief functions that deal with constructing and initialising the beam-column element
     */
    //@{
        // all element initialisers are element-specific.
    //@}
    /**
     * @name element operator overloads
     * @brief allows for overloading some operators for easier operations and for some common STL algorithms such as for sorting.
     * 
     */
    //@{
    /**
     * @brief overloads the less than operator to compare elements by their ID, allowing easy sorting of element STL containers via \ref std::sort.
     * 
     */
    virtual bool operator<(const BeamElementCommonInterface& other_elem) const
    { 
        return id < other_elem.id; 
    } 
    //@}
    /**
     * @name element property calculation functions
     * @brief functions that are used to calculate element properties such as stress, shape, length, and strain
     */
    //@{
        // all element property calculation functions are element-specific.
    //@}

    /**
     * @name stiffness matrix functions
     * @brief functions that deal with generating and evaluating the different stiffness matrices.
     */
    //@{
        /**
         * @brief runs all stiffness function calculations
         * 
         */
        virtual void calc_stiffnesses()
        {
            calc_mat_stiffness();
            calc_geom_stiffness();
            calc_tangent_stiffness();
            calc_elem_global_stiffness();
        }
    //@}

    /**
     * @name logging functions
     * @brief functions used for logging output to stream - mostly for debugging.
     */
    //@{
        /**
         * @brief prints the most important information of the element to the output stream.
         */
        virtual void print_info() {
            std::cout << "elem " << id << " of type " <<elem_type << " with " << ndofs << " dofs, and " << nnodes << " nodes:" << std::endl;
            for (auto node_i: nodes) {
                node_i->print_info();
            }
            std::cout << "it is also of length " << length << std::endl;
        }
         /**
         * @brief prints the internal state of the element.
         * 
         */
        virtual void print_element_state(bool print_nodal_disp = false, bool print_strains = false, 
                                            bool print_stresses = true, bool print_nodal_forces = false) 
        {
            if (print_nodal_disp) 
            {
                std::cout << "element " << id << " nodal displacements are:" << std::endl << local_d << std::endl;
            }
            if (print_strains) 
            {
                std::cout << "element " << id << " strains are:"  << std::endl << local_eps << std::endl;
            }
            if (print_stresses) 
            {
                std::cout << "element " << id << " stresses are:"  << std::endl << local_stresses << std::endl;
            }
            if (print_nodal_forces) 
            {
                std::cout << "element " << id << " nodal forces are:"  << std::endl << local_f << std::endl;
            }

        }
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
        virtual void get_U_from_nodes() 
        {
            std::array<real, 6> nodal_disp;
            int i = 0;
            // global_ele_U
            for (auto node: nodes)
            {
                nodal_disp = node->get_nodal_displacements();
                for (auto dof: nodal_disp)
                {
                    global_ele_U(i) = dof;
                    ++i;
                }
            }
        }

        /**
         * @brief Populates the resistance forces triplets removing any inactive freedoms.
         * @todo I seem to be doing the active_dofs thing too often. Perhaps the element should also have a set that contains its active dofs?
         */
        virtual void populate_resistance_force_triplets() {
            global_R_triplets.clear();
            // the 12x1 full resistance vector from local nodal forces vector f
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
                    force_value = element_global_resistance_forces(active_dof + total_nodal_ndofs_completed);
                    // since inactive nodes do not appear in R, we have to make sure to be careful about where we add our nodal forces.
                    // here, nz_i + nodal_dof_index simply starts at where the node freedoms start in the global index, and then
                    // iterates one by one. See how we ++ nodal_dof_index for each freedom we add, and how we restart from zero when
                    // we start work with the next node?
                    global_R_triplets.push_back(spnz(nz_i + nodal_dof_index, 0, force_value));
                    nodal_dof_index++;
                }
                //**< has to be 6 because each node has 6 dofs and our \ref element_global_resistance_forces also has 6 rows for each node!*
                total_nodal_ndofs_completed += 6;
            }
        }
        /**
         * @brief calculates the global stiffness contribution of the local element and populates global_stiffness_triplets
         * 
         * @details first, the freedoms are mapped to the right size getting the number of elements of the \ref elem_global_stiffness
         * After that, \ref stiffness_map is used to map where these contributions would go in the global stiffness
         * matrix. So, this function will populate \ref global_stiffness_triplets with sparse matrix notation
         * 
         */
        virtual void calc_global_stiffness_triplets() 
        {
            global_stiffness_triplets.clear();
            // we have the same number of contribution as stiffness components 
            // assuming all are non-zero, which is not correct but okay as it is "safe" although not very memory efficient.
            global_stiffness_triplets.reserve(elem_global_stiffness.rows() * elem_global_stiffness.cols());
            for (auto kmap: stiffness_map)
            {
                real val = elem_global_stiffness(kmap[0], kmap[1]);
                global_stiffness_triplets.push_back(spnz(kmap[2], kmap[3], val));
            }
        }

        /**
         * @brief populates \ref stiffness_map considering active and inactive DOFs for each node of the element
         * 
         * @details see function \ref calc_global_stiffness_triplets, and variables \ref stiffness_map, and \ref global_stiffness_triplets. 
         * 
         * @todo REALLY needs to be revisited. attempt to rewrite this function so it does the following:
         *  1. gets all the contribution without worrying about active or not
         *  2. if a contribution is inactive then that contribution is zeroed AND
         *  3. zeroed contributions are not added to \ref global_stiffness_triplets
         * 
         */
        virtual void map_stiffness()
        {
             // local to global stiffness map: <<local_row, local_col, global_row, global_col>, ...>
            stiffness_map.clear();
            int stiffness_size = 0;
            for (auto node: nodes) 
            {
                stiffness_size += std::size(node->get_active_dofs());
            }
            stiffness_size *= stiffness_size;
           
            stiffness_map.reserve(stiffness_size);
            int i = 0;
            for (auto node_i: nodes)
            {
                int j = 0;
                std::set<int> active_dofs_i = node_i->get_active_dofs();
                int nz_i_i = node_i->get_nz_i();
                for (auto node_j: nodes)
                {
                    
                    std::set<int> active_dofs_j = node_j->get_active_dofs();
                    
                    int nz_i_j = node_j->get_nz_i();
                    int dof_i_index = 0;
                    for (auto dof_i: active_dofs_i)
                    {
                        int dof_j_index = 0;
                        for (auto dof_j: active_dofs_j)
                        {
                            stiffness_map.push_back({6*i+dof_i, 6*j+dof_j, nz_i_i + dof_i_index, nz_i_j+dof_j_index});
                            ++dof_j_index;
                        }
                        ++dof_i_index;
                    }
                ++j;
                }
            ++i;
            }
        }
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
        virtual void set_global_U(vec global_U_vec) {global_ele_U = global_U_vec;}
        /**
         * @brief Set \ref local_d to some displacement vector.
         * 
         * @param new_disp the new displacement the \ref local_d would be replaced by.
         */
        virtual void set_d(vec new_disp) {local_d = new_disp;}
        
    //@}
    /**
     * @name getter functions
     * @brief functions that retrieve protected variables
     */
    //@{
        virtual int get_ndofs() const {return ndofs;}
        virtual int get_nnodes() const {return nnodes;}
        virtual std::string get_elem_type() const {return elem_type;}
        virtual unsigned get_id() const {return id;}

        virtual vec get_global_ele_U() const {return global_ele_U;}
        virtual vec get_local_d() const {return local_d;}
        virtual vec get_local_f() const {return local_f;}
        virtual vec get_element_resistance_forces() const {return element_global_resistance_forces;}   
        virtual std::vector<spnz> get_global_resistance_force_triplets() {return global_R_triplets;}        
        virtual vec get_eps() const {return local_eps;}
        virtual vec get_local_stresses() const {return local_stresses;}
        virtual mat get_local_constitutive_mat() const {return local_constitutive_mat;}
        virtual mat get_local_mat_stiffness() const {return local_mat_stiffness;}
        virtual mat get_local_geom_stiffness() const {return local_geom_stiffness;}
        virtual mat get_local_tangent_stiffness() const {return local_tangent_stiffness;}
        virtual mat get_elem_global_stiffness() const {return elem_global_stiffness;}
        virtual std::vector<spnz> get_global_stiffness_triplets() {return global_stiffness_triplets;}
        
        virtual mat get_N() const {return N;}
        virtual mat get_B() const {return B;}
        virtual mat get_T() {return orient.get_T();}
        virtual real get_L() {return orient.get_length();}

        virtual int const get_nth_node_id(int n) const 
        {
            if (n > nnodes - 1 || n < 0)
            {
                std::cout << "Error: Requested invalid node " << n << " from element " << id << std::endl;
                std::cout << "Element has " << nnodes << " nodes." << std::endl;
                std::exit(1);
            }
            return nodes[n]->get_id();
        }
    //@}
};

#endif