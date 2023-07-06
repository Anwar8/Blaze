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
        int nz_i = 0;
        std::set<int> connected_elements;
        std::set<int> active_dofs = {0, 1, 2, 3, 4, 5};
        std::set<int> inactive_dofs;
    public:
        Node();
        Node(real x_pos, real y_pos, real z_pos);
        Node(int i, coords xyz);
        void print_info();
        coords const get_coords() const;
        int const get_ndof()  {return ndof;}
        void add_connected_element(int element_id) {connected_elements.insert(element_id);}
        unsigned const get_id() const {return id;}
        
        // functions to work on DoFs vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
        std::set<int> const get_inactive_dofs() const {return inactive_dofs;}
        std::set<int> const get_active_dofs() const {return active_dofs;}
        void calc_ndof() {ndof = std::size(active_dofs);};
        bool valid_dof(int dof) {return (dof >= 0 && dof < 6);}
        void fix_dof(int dof);
        void free_dof(int dof);
        void set_nz_i(int i) {nz_i = i;}
        int get_nz_i() {return nz_i;}
        

        template <typename STLContainer>
        void fix_dofs(STLContainer dofs) {
            for (auto dof: dofs) {
                fix_dof(dof);
            }
            calc_ndof();
        }
        template <typename STLContainer>
        void free_dofs(STLContainer dofs) {
            for (auto dof: dofs) {
                free_dof(dof);
            }
            calc_ndof();
        }
        
        void fix_all_dofs();

        void free_all_dofs();
        void print_inactive_dofs();
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