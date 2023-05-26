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


};


#endif