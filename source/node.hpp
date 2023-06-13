#ifndef NODE_HPP
#define NODE_HPP
#include <set>
#include "maths_defaults.hpp"

class Node {
    private:
        unsigned id = 0;
        coords coordinates;
        real mass;
        int ndof = 6;
        std::set<int> connected_elements;
        std::set<int> inactive_dofs;
    public:
        Node();
        Node(real x_pos, real y_pos, real z_pos);
        Node(int i, coords xyz);
        void print_info();
        coords const get_coords() const;
        void add_connected_element(int element_id) {connected_elements.insert(element_id);}
        unsigned const get_id() const {return id;}
        
        // functions to work on DoFs vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
        void calc_ndof() {ndof = 6 - std::size(inactive_dofs);};
        bool valid_dof(int dof) {return (dof >= 0 && dof < 6);}
        void fix_dof(int dof);
        void free_dof(int dof);
        

        template <typename STLContainer>
        void fix_dofs(STLContainer dofs) {
            for (auto dof: dofs) {
                fix_dof(dof);
            }
        }
        template <typename STLContainer>
        void free_dofs(STLContainer dofs) {
            for (auto dof: dofs) {
                free_dof(dof);
            }
        }
        
        void fix_all_dofs() {inactive_dofs.insert({0, 1, 2, 3, 4, 5});}
        void free_all_dofs() {inactive_dofs.clear();}
        // functions to work on DoFs ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

        void  set_z(real z) { coordinates[2] = z;}
};

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