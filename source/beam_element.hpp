#ifndef BEAM_ELEMENT_HPP
#define BEAM_ELEMENT_HPP

#include <string>
#include <array>
#include "node.hpp"

class beam_element {
    private:
        unsigned id = 0;
        std::string const elem_type = "beam-column";
        int const dofs = 3;
        int const nnodes = 2;
        std::array<node, 2> nodes;

        double length = 0.0;
    public:
        beam_element();
        beam_element(std::array<node, 2> input_nodes);
        void print_info();
};
#endif