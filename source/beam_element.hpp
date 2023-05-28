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
        real I1 = 0.0004570000; // m^4
        real I2 = 0.0000235000; // m^4 
    public:
        real get_E() {return E;}
        real get_A() {return A;}
        real get_I1() {return I1;}
        real get_I2() {return I2;}
};

class shape_function {
    private:
        mat K = make_xd_mat(6,6);
        mat N = make_xd_mat(2,6);
        mat B = make_xd_mat(2,6);
    public:
        mat get_K() {return K;}
        mat get_N() {return N;}
        mat get_B() {return B;}
        void calc_N(real x, real L);
        void calc_B(real x, real L);
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
        void calc_eps();
        mat get_N() {return shape_func.get_N();}
        mat get_B() {return shape_func.get_B();}
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