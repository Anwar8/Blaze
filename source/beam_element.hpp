#ifndef BEAM_ELEMENT_HPP
#define BEAM_ELEMENT_HPP

#include <string>
#include <array>
#include "maths_defaults.hpp"
#include "node.hpp"

class beam_section {
    // 457 x 191 x 98
    private:
        real E = 2.06e11; // Pa
        real A = 0.0125; // m^2
        real I = 0.0004570000; // m^4
    public:
        real get_E() {return E;}
        real get_A() {return A;}
        real get_I() {return I;}
};

class shape_function {
    private:
        mat K = make_xd_mat(6,6);
        mat N = make_xd_mat(2,6);
        mat B = make_xd_mat(2,6);
    public:
        mat const get_K() const {return K;}
        mat const get_N() const {return N;}
        mat const get_B() const {return B;}
        void calc_N(real x, real L);
        void calc_B(real x, real L);
        void calc_K(real L, beam_section& sec);
};

class orientation {
    private:
        coords local_x;
        real length;
        mat T = make_xd_mat(6,6);
        real alpha = 0.0;
    public:
        void evaluate(std::array<node, 2> const & nodes, coords const & origin_x)
        {
            calc_length_local_x(nodes);
            calc_alpha(origin_x, nodes[0].get_coords());
            calc_T();
        }
        void calc_length_local_x(std::array<node, 2> const &  nodes) {
            local_x = (nodes[1].get_coords() - nodes[0].get_coords());
            length = local_x.norm();
            local_x /= length;
        }
        void calc_alpha(coords const& origin_x, coords const node_i_coords) {
            // calcualte the angle between the local and global x axes
           alpha = std::acos(origin_x.dot(local_x));
        }
        void calc_T() {
            T(0,0) = std::cos(alpha);
            T(1,0) = -std::sin(alpha);
            T(0,1) = std::sin(alpha);
            T(1,1) = std::cos(alpha);
            T(2,2) = 1;
            T(3,3) = std::cos(alpha);
            T(4,3) = -std::sin(alpha);
            T(3,4) = std::sin(alpha);
            T(4,4) = std::cos(alpha);
            T(5,5) = 1;
        }
        mat get_T() {return T;}
};

class beam_element {
    private:
        unsigned id = 0;
        std::string const elem_type = "beam-column";
        int const dofs = 3;
        int const nnodes = 2;
        std::array<node, 2> nodes;
        
        beam_section section;
        shape_function shape_func;
        orientation orient;
        real length = 0.0;
        vec local_d = make_xd_vec(6);
        vec local_f = make_xd_vec(6);
        vec local_eps = make_xd_vec(2);

    public:
        beam_element();
        beam_element(std::array<node, 2> input_nodes);
        void print_info();
        void calc_length();
        void calc_N(real x);
        void calc_B(real x);
        void calc_K();
        void calc_T(coords origin_x);
        void calc_eps();
        mat get_N() {return shape_func.get_N();}
        mat get_B() {return shape_func.get_B();}
        mat get_K() {return shape_func.get_K();}
        mat get_T() {return orient.get_T();}
        vec get_eps() {return local_eps;}
        vec get_d() {return local_d;}

        void set_d(vec new_disp) {local_d = new_disp;}

    beam_element(const beam_element& other) {
        id = other.id;
        nodes = other.nodes;
        section = other.section;
        shape_func = other.shape_func;
        length = other.length;
        local_d = other.local_d;
        local_f = other.local_f;
        local_eps = other.local_eps;
    }
    // Copy assignment - i need to learn a bit more...
    beam_element& operator=(const beam_element& other) {
        if (this != &other) {
            beam_element temp(other); // use copy constructor to create temporary object

            std::swap(id, temp.id);
            std::swap(nodes, temp.nodes);
            std::swap(section, temp.section);
            std::swap(shape_func, temp.shape_func);
            std::swap(length, temp.length);
            std::swap(local_d, temp.local_d);
            std::swap(local_f, temp.local_f);
            std::swap(local_eps, temp.local_eps);
        }
        return *this;
    }

};


#endif