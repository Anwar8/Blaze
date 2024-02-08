#ifndef NODE_HPP
#define NODE_HPP
#include <set>
#include "maths_defaults.hpp"


/**
 * @brief node data and functions including functionality to activate and deactivate DOFs.
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
        int ndof = 6; /**< the number of DOFs - should be 6 unless some are deactivated.*/
        int nz_i = 0; /**< FORGOT - used in assembly.*/
        std::set<int> connected_elements; /**< of element ids that are connected to this node; expected to be useful for element and node deletion.*/
        std::set<int> active_dofs = {0, 1, 2, 3, 4, 5}; /**< set of active DOFs; all of them at first, then if deactivated moved to inactive_dofs.*/
        std::set<int> inactive_dofs; /**< a std set of active DOFs; none at first, then if any are deactivated then they are moved from active_dofs*/
    public:
        /**
         * @brief Construct a new Node object with 0 mass and 0 across coordinates
         * 
         */
        Node();

        /**
         * @brief Construct a new Node object with 0 mass and specified x, y, and z coordinates
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
         * @name setters_getters
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
         * @attention \ref nz_i is still a mystery. 
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
         * @brief checks if a given DOF is valid within the problem domain which has 6 DOFs by default.
         * 
         * @param dof  DOF to check.
         * @return true if DOF is valid and is between 0 and 5.
         * @return false if DOF is not a valid DOF (not between 0 and 5).
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
         * @param dof degree of freedom to activate (free).
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
};

/**
 * @brief a class to contain the global coordinate system using unit vectors.
 * 
 * @details The global coordinate system can allow members to transform
 * themselves with respect to it using their /ref BasicOrientation object.
 * 
 */
class GlobalCoords {
    private:
        coords centroid = {0.0, 0.0 , 0.0};
        coords unit_x = {1.0, 0.0, 0.0};
        coords unit_y = {0.0, 1.0, 0.0};
        coords unit_z = {0.0, 0.0, 1.0};
    public: 
        coords get_centroid() {return centroid;}
        coords get_unit_x() {return unit_x;}
        coords get_unit_y() {return unit_y;}
        coords get_unit_z() {return unit_z;}
};

#endif