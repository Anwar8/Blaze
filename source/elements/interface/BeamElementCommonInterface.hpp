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
template <typename BeamSectionClass>
class BeamElementCommonInterface : public BeamElementBaseClass<BeamSectionClass> {
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
        return this->id < other_elem.id; 
    } 
    //@}
    /**
     * @name element property calculation functions
     * @brief functions that are used to calculate element properties such as stress, shape, length, and strain
     */
    //@{
        /**
         * @brief updates the starting state of the section after solution convergence. Please see Bhatti's Advanced Topics in Finite Element Analysis of Structures for more on this. Does nothing for all Elastic elements.
         */
        virtual void update_section_starting_state() override {};
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
        virtual void calc_stiffnesses() override
        {
            this->calc_mat_stiffness();
            this->calc_geom_stiffness();
            this->calc_tangent_stiffness();
            this->calc_external_geom_stiffness();
            this->calc_elem_global_stiffness();
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
        virtual void print_info() override {
            std::cout << "elem " << this->id << " of type " <<this->elem_type << " with " << this->ndofs << " dofs, and " << this->nnodes << " nodes:" << std::endl;
            for (auto& node_i: this->nodes) {
                node_i->print_info();
            }
            std::cout << "it is also of length " << this->length << std::endl;
        }
         /**
         * @brief prints the internal state of the element.
         * 
         */
        virtual void print_element_state(bool print_nodal_disp = false, bool print_strains = false, 
                                            bool print_stresses = true, bool print_nodal_forces = false) override
        {
            if (print_nodal_disp) 
            {
                std::cout << "element " << this->id << " nodal displacements are:" << std::endl << this->local_d << std::endl;
            }
            if (print_strains) 
            {
                std::cout << "element " << this->id << " strains are:"  << std::endl << this->local_eps[0] << std::endl;
            }
            if (print_stresses) 
            {
                std::cout << "element " << this->id << " stresses are:"  << std::endl << this->local_stresses[0] << std::endl;
            }
            if (print_nodal_forces) 
            {
                std::cout << "element " << this->id << " nodal forces are:"  << std::endl << this->local_f << std::endl;
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
        virtual void get_U_from_nodes() override
        {
            std::array<real, 6> nodal_disp;
            int i = 0;
            // global_ele_U
            for (auto& node: this->nodes)
            {
                nodal_disp = node->get_nodal_displacements();
                for (auto& dof: nodal_disp)
                {
                    this->global_ele_U(i) = dof;
                    ++i;
                }
            }
        }

        /**
         * @brief Populates the resistance forces triplets removing any inactive freedoms.
         * @todo I seem to be doing the active_dofs thing too often. Perhaps the element should also have a set that contains its active dofs?
         */
        virtual void populate_resistance_force_triplets() override {
            this->global_R_triplets.clear();
            // the 12x1 full resistance vector from local nodal forces vector f
            std::set<int> node_active_dofs;
            int nz_i = 0;
            real force_value;
            int total_nodal_ndofs_completed = 0; // each node we finish with, we add 6 to this. 
            // This means we have to move to the next set of values corresponding to the next 
            // node in the full resistance vector.         
            for (auto& node: this->nodes)
            {
                int nodal_dof_index = 0;
                node_active_dofs = node->get_active_dofs();
                nz_i = node->get_nz_i();
                for (auto& active_dof: node_active_dofs)
                {
                    force_value = this->element_global_resistance_forces(active_dof + total_nodal_ndofs_completed);
                    // since inactive nodes do not appear in R, we have to make sure to be careful about where we add our nodal forces.
                    // here, nz_i + nodal_dof_index simply starts at where the node freedoms start in the global index, and then
                    // iterates one by one. See how we ++ nodal_dof_index for each freedom we add, and how we restart from zero when
                    // we start work with the next node?
                    this->global_R_triplets.push_back(spnz(nz_i + nodal_dof_index, 0, force_value));
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
        virtual void calc_global_stiffness_triplets() override
        {
            this->global_stiffness_triplets.clear();
            // we have the same number of contribution as stiffness components 
            // assuming all are non-zero, which is not correct but okay as it is "safe" although not very memory efficient.
            this->global_stiffness_triplets.reserve(this->stiffness_map.size());
            for (auto& kmap: this->stiffness_map)
            {
                real val = this->elem_global_stiffness(kmap[0], kmap[1]);
                this->global_stiffness_triplets.push_back(spnz(kmap[2], kmap[3], val));
            }
        }
        
        /**
         * @brief populates \ref stiffness_map considering active and inactive DOFs for each node of the element. Uses the position vectors of the element to do this. Does not collect row-contributions from interface nodes i.e. those not on their parent rank.
         * 
         * @details see function \ref calc_global_stiffness_triplets, and variables \ref stiffness_map, and \ref global_stiffness_triplets. 
         * 
         */
        virtual void map_stiffness() override
        {
            std::vector<int> local_position_vector_rows;
            std::vector<int> local_position_vector_columns;
            std::vector<int> position_vector_rows;
            std::vector<int> position_vector_columns;

            std::vector<int> previous_nodes_ndofs;
            int element_ndofs = 0;
            for (auto& node: this->nodes) 
            {
                previous_nodes_ndofs.push_back(element_ndofs);
                element_ndofs += node->get_ndof();
            }
            
            // set the elements of the position_vector_for_rows
            for (auto& node: this->nodes) 
            {
                std::vector<int> node_dofs_numbers;
                if (node->is_on_parent_rank())
                {
                    node_dofs_numbers = node->get_dofs_numbers();
                    position_vector_rows.insert(position_vector_rows.end(), node_dofs_numbers.begin(), node_dofs_numbers.end());
                }
            }
            // set the elements of the position_vector_for_columns
            for (auto& node: this->nodes) 
            {
                std::vector<int> node_dofs_numbers;
                node_dofs_numbers = node->get_dofs_numbers();
                position_vector_columns.insert(position_vector_columns.end(), node_dofs_numbers.begin(), node_dofs_numbers.end());
            }
            
            // set the elements of the local_position_vector_for_rows
            for (int node_i = 0; node_i < this->nnodes; ++node_i)
            {
                std::vector<int> node_active_dofs;
                std::set<int> node_active_dofs_set;
                if ( this->nodes[node_i]->is_on_parent_rank())
                {
                    node_active_dofs_set = this->nodes[node_i]->get_active_dofs();
                    node_active_dofs = std::vector<int>(node_active_dofs_set.begin(), node_active_dofs_set.end());
                    std::for_each(node_active_dofs.begin(), node_active_dofs.end(), [=](int& active_dof) { active_dof+=node_i*6;}); 
                    local_position_vector_rows.insert(local_position_vector_rows.end(), node_active_dofs.begin(), node_active_dofs.end());
                }
            }

            // set the elements of the local_position_vector_for_columns
            for (int node_i = 0; node_i < this->nnodes; ++node_i)
            {
                std::vector<int> node_active_dofs;
                std::set<int> node_active_dofs_set;
                node_active_dofs_set = this->nodes[node_i]->get_active_dofs();
                node_active_dofs = std::vector<int>(node_active_dofs_set.begin(), node_active_dofs_set.end());
                std::for_each(node_active_dofs.begin(), node_active_dofs.end(), [=](int& active_dof) { active_dof+=node_i*6;}); 
                local_position_vector_columns.insert(local_position_vector_columns.end(), node_active_dofs.begin(), node_active_dofs.end());
            }

            this->stiffness_map.clear();
            int stiffness_size = local_position_vector_rows.size() * local_position_vector_columns.size();
            this->stiffness_map.reserve(stiffness_size);

            for (int row_i  = 0; row_i < local_position_vector_rows.size(); ++row_i)
            {
                for (int column_i = 0; column_i < local_position_vector_columns.size(); ++column_i)
                {
                    this->stiffness_map.push_back({local_position_vector_rows[row_i], local_position_vector_columns[column_i], position_vector_rows[row_i], position_vector_columns[column_i]});
                }
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
        virtual void set_global_U(vec global_U_vec) override {this->global_ele_U = global_U_vec;}
        /**
         * @brief Set \ref local_d to some displacement vector.
         * 
         * @param new_disp the new displacement the \ref local_d would be replaced by.
         */
        virtual void set_d(vec new_disp) override {this->local_d = new_disp;}
        
    //@}
    /**
     * @name getter functions
     * @brief functions that retrieve protected variables
     */
    //@{
        virtual int get_ndofs() const override {return this->ndofs;}
        virtual int get_nnodes() const override {return this->nnodes;}
        virtual std::string get_elem_type() const override {return this->elem_type;}
        virtual unsigned get_id() const override {return this->id;}

        virtual vec get_global_ele_U() const override {return this->global_ele_U;}
        virtual vec get_local_d() const override {return this->local_d;}
        virtual vec get_local_f() const override {return this->local_f;}
        virtual vec get_element_resistance_forces() const override {return this->element_global_resistance_forces;}   
        virtual std::vector<spnz> get_global_resistance_force_triplets() override {return this->global_R_triplets;}  
        /**
         * @brief inserts the contents of \ref element_global_resistance_forces into the end of global resistance triplets vector \ref Assembler::R_global_triplets. 
         * Used to reduce copying during assembly, still basically a getter function.
         * 
         * @param R_global_triplets The container for the triplets that are used for assembling the global resistance vector \f$ \boldsymbol{R}\f$ from element contributions.
         */ 
        virtual void insert_global_resistance_force_triplets(std::vector<spnz>& global_resistance_triplets_vector) override 
        {
            global_resistance_triplets_vector.insert(global_resistance_triplets_vector.end(), this->global_R_triplets.begin(), this->global_R_triplets.end());            
        }      
        virtual vec get_eps() const override {return this->local_eps[0];}
        virtual vec get_local_stresses() const override {return this->local_stresses[0];}
        virtual mat get_local_constitutive_mat() const override {return this->local_constitutive_mat[0];}
        virtual mat get_local_mat_stiffness() const override {return this->local_mat_stiffness;}
        virtual mat get_local_geom_stiffness() const override {return this->local_geom_stiffness;}
        virtual mat get_local_tangent_stiffness() const override {return this->local_tangent_stiffness;}
        virtual mat get_elem_global_stiffness() const override {return this->elem_global_stiffness;}
        virtual std::vector<spnz> get_global_stiffness_triplets() override {return this->global_stiffness_triplets;}
        /**
         * @brief inserts the contents of \ref global_stiffness_triplets into the end of global_triplets_vector. Used to reduce copying during assembly, still basically a getter function.
         * 
         * @param global_triplets_vector The container for the triplets that are used for assembling the global stiffness matrix \f$ \boldsymbol{K}\f$ from element contributions.
         */
        virtual void insert_global_stiffness_triplets(std::vector<spnz>& global_triplets_vector) override {
            global_triplets_vector.insert(global_triplets_vector.end(), this->global_stiffness_triplets.begin(), this->global_stiffness_triplets.end());
        }
        virtual mat get_N() const override {return this->N[0];}
        virtual mat get_B() const override {return this->B[0];}
        virtual mat get_T() override {return this->orient.get_T();}
        virtual real get_L() override {return this->orient.get_length();}

        virtual int const get_nth_node_id(int n) const override 
        {
            if (n > this->nnodes - 1 || n < 0)
            {
                std::cout << "Error: Requested invalid node " << n << " from element " << this->id << std::endl;
                std::cout << "Element has " << this->nnodes << " nodes." << std::endl;
                std::exit(1);
            }
            return this->nodes[n]->get_id();
        }
    //@}
};

#endif