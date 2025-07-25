/**
 * @file node.hpp
 * @brief node object and degree of freedom controls
 * 
 */

#ifndef NODE_HPP
#define NODE_HPP
#include <set>
#include "basic_utilities.hpp"
#include "maths_defaults.hpp"
#include "blaze_config.hpp"


/**
 * @brief node data and functions including functionality to activate and deactivate DoFs.
 * 
 * @todo Add functionality and interface to apply mechanical loading on nodes.
 * 
 */
class Node {

    private:
        unsigned id = 0; /**< unique id for node object - can be changed during renumbering.*/
        unsigned record_id = 0; /**< unique id for node object - cannot be changed from original value given in mesh.*/
        coords coordinates; /**< x, y, and z coordinates of the node.*/
        /**
         * @brief nodal mass.
         * 
         * @todo nodal displacements should be Eigen3 vector types not std array types.
         */
        real mass;
        int ndof = 6; /**< the number of DoFs - should be 6 unless some are deactivated.*/
        int nz_i = 0; /**< Corresponds to global location of node and its DoFs and load, considering deactivated DoFs.*/
        bool on_parent_rank = true; /**< Is this node object currently living on its parent rank?*/
        int parent_rank = 0; /**< The rank that owns this node.*/
        std::set<int> connected_elements; /**< set of element ids that are connected to this node; expected to be useful for element and node deletion.*/
        std::set<int> active_dofs = {0, 1, 2, 3, 4, 5}; /**< set of active DOFs; all of them at first, then if deactivated moved to inactive_dofs.*/
        std::set<int> inactive_dofs; /**< a std set of active DoFs; none at first, then if any are deactivated then they are moved from active_dofs.*/
        std::vector<int> dofs_numbers = {0, 1, 2, 3, 4, 5}; /**< a std vector of DoF numbers where each number is equal to increments on nz_i; equal to \ref active_dofs at first, then updated when \ref nz_i is known.*/
        
        std::set<int> loaded_dofs; /**< a std set of loaded DoFs; none at first, then those loaded are added.*/
        std::array<real, 6> nodal_loads = {0., 0., 0., 0., 0., 0.}; /**< a std array containing 6 slots to be filled with nodal loads corresponding to dofs; initialised to zero.*/
        std::vector<spnz> global_nodal_loads_triplets; /**< the global contributions of the element to the global stiffness - made as sparse matrix contributions that would be gathered to create the global sparse matrix*/
        
        std::array<real, 6> nodal_displacements = {0., 0., 0., 0., 0., 0.}; /**< a std array containing 6 slots to be filled with nodal displacements corresponding to dofs; initialised to zero.*/
    public:
        /**
         * @brief Construct a new Node object with 0 mass and 0 across coordinates.
         * 
         */
        Node() : coordinates(0.0 , 0.0, 0.0), mass(0.0) {}

        /**
         * @brief Construct a new Node object with 0 mass and specified x, y, and z coordinates.
         * 
         * @attention uses the default id of 0; pretty pointless but useful during development and testing.
         * 
         * @param x_pos x coordinate.
         * @param y_pos y coordinate.
         * @param z_pos z coordinate.
         */
        Node(real x_pos, real y_pos, real z_pos) : coordinates(x_pos, y_pos, z_pos), mass(0.0) {}

        /**
         * @brief Construct a new Node object with 0 mass, and an id and coordinates.
         * 
         * @param i id number.
         * @param xyz nodal x, y, and z coordinates.
         */
        Node(int i, coords xyz) : id(i), record_id(i), coordinates(xyz), mass(0.0) {}

        /**
         * @brief overloads the less than operator to compare nodes by their node ID, allowing easy sorting of node STL containers via \ref std::sort.
         * 
         */
        bool operator<(const Node& other_node) const
        { 
            return id < other_node.id; 
        } 

        void print_info()
        {
            std::cout << "Node " << id << ": xyz = (" << coordinates[0] << ", " << coordinates[1] << ", " << coordinates[2] <<  "), and mass = " << mass << std::endl;
            std::cout << "There are " << std::size(connected_elements) << " connected elements. They are: ";
            print_container<std::set<int>>(connected_elements);
            std::cout << "Node has following loads:" << std::endl;
            print_container(nodal_loads);
            std::cout << "Node has following displacement:" << std::endl;
            print_container(nodal_displacements);
            print_inactive_dofs();
        }
        /**
         * @name SettersGetters
         * @brief functions to set and get private class members.
         */
        //@{
        /**
         * @brief Get the coordinates of the node.
         * 
         * @todo inline this function so that it directly returns the coordinates.
         * 
         * @return coords const coordinates of the node.
         */
        coords const get_coords() const {return coordinates;}
        int const get_ndof()  {return ndof;}
        void add_connected_element(int element_id) {connected_elements.insert(element_id);}
        unsigned const get_id() const {return id;}
        unsigned const get_record_id() const {return record_id;}
        int const get_num_connected_elements() {return connected_elements.size();}
        /**
         * @brief Get the number of rows of the stiffness matrix to which this node contributes. The calculation assumes that the self-contributions of the node to the stiffness matrix sum into the same spot for all elements connected to the node after the first one, leaving only the stiffness components that are connected to other nodes. The reason we have 2*ndof is because each element has 2 nodes.
         * @warning only works assuming each element has two nodes!
         * @return int const 
         */
        int const get_num_row_contributions() 
        {
            return ndof + get_num_connected_elements()*ndof;
        }
        /**
         * @brief Set the nz_i to a value, and calls \ref update_dofs_numbers.
         * 
         * @attention \ref nz_i global index of node considering activated and deactivated DoFs.
         * 
         * @param i the value to set nz_i to; can be any integer.
         */
        void set_nz_i(int i) {
            nz_i = i;
            update_dofs_numbers();
        }
        /**
         * @brief Ipdates all members of \ref dofs_numbers due to change in \ref nz_i. 
         * 
         */
        void update_dofs_numbers()
        {
            dofs_numbers.clear();
            for (int i = 0; i < active_dofs.size(); ++i)
            {
                dofs_numbers.push_back(nz_i + i);
            }
        }

        std::vector<int> get_dofs_numbers()
        {
            return dofs_numbers;
        }
        /**
         * @brief Increments the nz_i by a value, and calls \ref update_dofs_numbers.
         *
         * @param i the value to increment nz_i by; can be any integer.
         */
        void increment_nz_i(int i) {
            nz_i += i;
            update_dofs_numbers();
        }
        int get_nz_i() {return nz_i;}
        void  set_z(real z) { coordinates[2] = z;}

        void set_id(unsigned new_id) { id = new_id; }
        void increment_id(unsigned id_increment) { id += id_increment; }
        void set_parent_rank(int parent_rank, int calling_rank) 
        {
            parent_rank = parent_rank;
            on_parent_rank = (parent_rank == calling_rank);
        }
        int get_parent_rank() const {return parent_rank;}
        bool is_on_parent_rank() const {return on_parent_rank;}
        //@}

        /**
         * @name nodal_dof_functions
         * @brief functions that operate on nodal DoFs to activate and deactivate them.
         */
        //@{
        std::set<int> const get_inactive_dofs() const {return inactive_dofs;}
        std::set<int> const get_active_dofs() const {return active_dofs;}
        /**
         * @brief calculates the number of DOFs by checking size of \ref active_dofs.
         * 
         */
        void calc_ndof() {ndof = std::size(active_dofs);};

        /**
         * @brief checks if a given DoF is valid within the problem domain which has 6 DoFs by default.
         * 
         * @param dof  DoF to check.
         * @return true if DoF is valid and is between 0 and 5.
         * @return false if DoF is not a valid DoF (not between 0 and 5).
         */
        bool valid_dof(int dof) {return (dof >= 0 && dof < 6);}
        /**
         * @brief deactivates a DoF.
         * 
         * @details inserts DoF into \ref inactive_dofs, and erases it from \ref active_dofs,
         * then caculates the number of active DoFs \ref ndof via \ref calc_ndof.
         * 
         * @param dof DoF to fix.
         */
        void fix_dof(int dof)
        {
            if (valid_dof(dof))
            {
                inactive_dofs.insert(dof);
                active_dofs.erase(dof);
                calc_ndof();
            } else {
                std::cout << "ERROR: Cannot fix DoF " << dof << ". Only DoFs 0 through 5 allowed." << std::endl;
                std::exit(1);
            }
        }

        /**
         * @brief activates a DoF.
         * 
         * @details Erases DoF from \ref inactive_dofs and inserts it into \ref active_dofs,
         * then caculates the number of active DoFs \ref ndof via \ref calc_ndof.
         * 
         * @attention does not care if the DoF is already active. If so, the usage 
         * of std set would take care of not having repeats in \ref inactive_dofs 
         * and \ref active_dofs.
         * 
         * @param dof degree of freedom to activate (to free).
         */
        void free_dof(int dof) {
            if (valid_dof(dof))
            {
                inactive_dofs.erase(dof);
                active_dofs.insert(dof);
                calc_ndof();
            } else {
                std::cout << "ERROR: Cannot free DoF " << dof << ". Only DoFs 0 through 5 allowed." << std::endl;
                std::exit(1);
            }
        }

        /**
         * @brief templated function to fix multiple DoFs using \ref fix_dof.
         * 
         * @tparam STLContainer any stl-compatible container with appropriate built-in iterators.
         * @param dofs the dofs to fix. Each must be a valid DoF between 0 and 5.
         */
        template <typename STLContainer>
        void fix_dofs(STLContainer dofs) {
            for (auto& dof: dofs) {
                fix_dof(dof);
            }
        }
        
        /**
         * @brief templated function to free multiple DoFs using \ref free_dof.
         * 
         * @tparam STLContainer any stl-compatible container with appropriate built-in iterators.
         * @param dofs the dofs to free. Each must be a valid DoF between 0 and 5.
         */
        template <typename STLContainer>
        void free_dofs(STLContainer dofs) {
            for (auto& dof: dofs) {
                free_dof(dof);
            }
        }
        
        void fix_all_dofs()
        {
            inactive_dofs.insert({0, 1, 2, 3, 4, 5});
            active_dofs.clear();
            calc_ndof();
        }
        void free_all_dofs()
        {
            inactive_dofs.clear();
            active_dofs.insert({0, 1, 2, 3, 4, 5});
            calc_ndof();
        }

        void print_inactive_dofs() 
        {
            std::cout << "Node " << id << " has " << std::size(inactive_dofs) << " inactive DoFs: ";
            print_container(inactive_dofs);
        }

        //@}
        /**
         * @name nodal_load_functions
         * @brief functions that operate on nodal DoF loads.
         */
        //@{      
        /**
         * @brief adds nodal load to \ref nodal_loads.
         * 
         * @param nodal_load nodal load to be added.
         * @param dof the DoF to which the nodal load will be added.
         */
        void add_nodal_load(real nodal_load, int dof) 
        {
            if (valid_dof(dof))
            {
                nodal_loads[dof] = nodal_load;
                loaded_dofs.insert(dof);

                if (VERBOSE)
                {
                std::cout << "node " << id << " loaded dofs are:" << std::endl;
                print_container(loaded_dofs);
                std::cout << "node " << id << " loads are:" << std::endl;
                print_container(nodal_loads);
                }
            } else {
                std::cout << "ERROR: Cannot add load to DoF " << dof << ". Only DoFs 0 through 5 can be loaded." << std::endl;
                std::exit(1);
            }
        }
        /**
         * @brief increments the nodal at DoF dof by dP.
         * 
         * @param nodal_load nodal load increment to be added to node loads.
         * @param dof the DoF to which the nodal load will be added.
         */
        void increment_nodal_load(real dP, int dof) 
        {
            if (valid_dof(dof))
            {
                if (loaded_dofs.count(dof)) {
                if (VERBOSE_NLB)
                {
                    std::cout << "Incrementing load at DoF " << dof << " of node " << id << " by " << dP << "." << std::endl;
                }
                nodal_loads[dof] += dP;
                } else {
                std::cout << "ERROR: Cannot increment load at DoF " << dof << "as this DoF is not already loaded." << std::endl;
                std::exit(1);
                }
            } else {
                std::cout << "ERROR: Cannot increment load to DoF " << dof << ". Only DoFs 0 through 5 can be loaded." << std::endl;
                std::exit(1);
            }
        }

        /**
         * @brief checks if the nodal loads are applied to inactive DoFs, and prints a warning if so.
         * 
         */
        void check_loads() 
        {
            for (auto& dof : loaded_dofs)
            {
                if (inactive_dofs.count(dof))
                {
                    std::cout << "WARNING: node "<< id << " DoF " << dof << " is inactive. A nodal load was added but will not be applied." << std::endl;
                }
            }
        }
        void clear_nodal_loads() 
        {
            nodal_loads = {0., 0., 0., 0., 0., 0.};
            loaded_dofs.clear();
        }
        /**
         * @brief converts the \ref nodal_loads array into a std vector of triplets to be collected by the assembler.
         * 
         * @warning requires C++20 or won't compile due to the use of the container.contains function introduced in the C++20 standard.
         */
        void compute_global_load_triplets()
        {
            global_nodal_loads_triplets.clear();
            int dof_index = 0;
            for (auto& active_dof: active_dofs) {
            if (VERBOSE)
            {
                std::cout << "node " << id << " checking active dof: " << active_dof << " with index " << dof_index << std::endl;
            }
            if (loaded_dofs.count(active_dof)) {
                if (VERBOSE)
                {
                std::cout << "pushing triplet val " << nodal_loads[active_dof] << " to P vector index " << nz_i + dof_index << std::endl;
                }
                global_nodal_loads_triplets.push_back(spnz(nz_i + dof_index, 0, nodal_loads[active_dof]));
            }
            dof_index++;
            }
        }
        /**
         * @brief returns the \ref global_nodal_loads_triplets vector.
         * 
         */
        std::vector<spnz> get_load_triplets() {return global_nodal_loads_triplets;}

        /**
         * @brief inserts the contents of \ref global_nodal_loads_triplets into the end of global_load_triplets_vector used to construct \f$\boldsymbol{P}\f$. Used to reduce copying during assembly, still basically a getter function.
         * 
         * @param global_load_triplets_vector The container for the triplets that are used for assembling the global load vector \f$ \boldsymbol{P}\f$ from nodal loads.
         */
        void insert_load_triplets(std::vector<spnz>& global_load_triplets_vector) {
            global_load_triplets_vector.insert(global_load_triplets_vector.end(), this->global_nodal_loads_triplets.begin(), this->global_nodal_loads_triplets.end());   
        }
        //@}
        
        /**
         * @name nodal_displacement_functions
         * @brief functions that operate on nodal DoF displacements.
         */
        //@{      
        /**
         * @brief sets the value of the nodal displacements.
         * 
         * @details sets a given DoF to a particular displacement value in the container \ref nodal_displacements.
         * 
         * @param dof degree of freedom to set the displacement value to.
         * @param disp the displacement value to which to set the DoF to.
         */
        void set_nodal_displacement(int dof, real disp)
        {
            if (valid_dof(dof)) {
                nodal_displacements[dof] = disp;
            }
        }
        /**
         * @brief Get the nodal displacements array.
         * 
         * @return std::array<real,6>  nodal_displacements array.
         */
        std::array<real,6> get_nodal_displacements() {return nodal_displacements;}
        /**
         * @brief Get the nodal displacements for a given dof.
         * 
         * @param dof the degree of freedom for which to return the displacement.
         * @return the nodal displacement corresponding to the dof requested.
         */
        real get_nodal_displacement(int dof) {return nodal_displacements[dof];}

        /**
         * @brief Get the loaded dofs set. Used for testing.
         * 
         * @return std::set<int> loaded_dofs.
         */
        std::set<int> get_loaded_dofs() {return loaded_dofs;}
        /**
         * @brief Get the nodal_loads std::array. Used for testing.
         * 
         * @return std::array<real, 6> nodal_loads.
         */
        std::array<real, 6> get_loads() {return nodal_loads;}
        //@}
};

#endif