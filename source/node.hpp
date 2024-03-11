/**
 * @file node.hpp
 * @brief node object and degree of freedom controls
 * 
 */

#ifndef NODE_HPP
#define NODE_HPP
#include <set>
#include "maths_defaults.hpp"


/**
 * @brief node data and functions including functionality to activate and deactivate DoFs.
 * 
 * @todo Add functionality and interface to apply mechanical loading on nodes.
 * 
 */
class Node {

    private:
        unsigned id = 0; /**< unique id for node object.*/
        coords coordinates; /**< x, y, and z coordinates of the node.*/
        /**
         * @brief nodal mass.
         * 
         * @todo should the nodal mass have x, y, and z components?
         */
        real mass;
        int ndof = 6; /**< the number of DoFs - should be 6 unless some are deactivated.*/
        int nz_i = 0; /**< Corresponds to global location of node and its DoFs and load, considering deactivated DoFs.*/
        std::set<int> connected_elements; /**< of element ids that are connected to this node; expected to be useful for element and node deletion.*/
        std::set<int> active_dofs = {0, 1, 2, 3, 4, 5}; /**< set of active DOFs; all of them at first, then if deactivated moved to inactive_dofs.*/
        std::set<int> inactive_dofs; /**< a std set of active DoFs; none at first, then if any are deactivated then they are moved from active_dofs.*/
        
        std::set<int> loaded_dofs; /**< a std set of loaded DoFs; none at first, then those loaded are added.*/
        std::array<real, 6> nodal_loads = {0., 0., 0., 0., 0., 0.}; /**< a std array containing 6 slots to be filled with nodal loads corresponding to dofs; initialised to zero.*/
        std::vector<spnz> global_nodal_loads; /**< the global contributions of the element to the global stiffness - made as sparse matrix contributions that would be gatehred to create the global sparse matrix*/
        
    public:
        /**
         * @brief Construct a new Node object with 0 mass and 0 across coordinates.
         * 
         */
        Node();

        /**
         * @brief Construct a new Node object with 0 mass and specified x, y, and z coordinates.
         * 
         * @attention uses the default id of 0; pretty pointless but useful during development and testing.
         * 
         * @param x_pos x coordinate.
         * @param y_pos y coordinate.
         * @param z_pos z coordinate.
         */
        Node(real x_pos, real y_pos, real z_pos);

        /**
         * @brief Construct a new Node object with 0 mass, and an id and coordinates.
         * 
         * @param i id number.
         * @param xyz nodal x, y, and z coordinates.
         */
        Node(int i, coords xyz);

        void print_info();

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
        coords const get_coords() const;
        int const get_ndof()  {return ndof;}
        void add_connected_element(int element_id) {connected_elements.insert(element_id);}
        unsigned const get_id() const {return id;}
        /**
         * @brief Set the nz_i to a value.
         * 
         * @attention \ref nz_i global index of node considering activated and deactivated DoFs.
         * 
         * @param i the value to set nz_i to; can be any integer.
         */
        void set_nz_i(int i) {nz_i = i;}
        int get_nz_i() {return nz_i;}
        void  set_z(real z) { coordinates[2] = z;}
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
        void fix_dof(int dof);

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
        void free_dof(int dof);

        /**
         * @brief templated function to fix multiple DoFs using \ref fix_dof.
         * 
         * @tparam STLContainer any stl-compatible container with appropriate built-in iterators.
         * @param dofs the dofs to fix. Each must be a valid DoF between 0 and 5.
         */
        template <typename STLContainer>
        void fix_dofs(STLContainer dofs) {
            for (auto dof: dofs) {
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
            for (auto dof: dofs) {
                free_dof(dof);
            }
        }
        
        void fix_all_dofs();
        void free_all_dofs();
        void print_inactive_dofs();

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
         * @param dof the dof to which the nodal load will be added.
         */
        void add_nodal_load(real nodal_load, int dof);
        /**
         * @brief converts the \ref nodal_loads array into a std vector of triplets to be collected by the assembler.
         * 
         * @warning requires C++20 or won't compile due to the use of the container.contains function introduced in the C++20 standard.
         */
        void compute_load_triplets();

};

#endif