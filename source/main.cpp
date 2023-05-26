#include <array>
#include "node.hpp"
#include "beam_element.hpp"

int main () {
    std::array<node, 2> in_nodes;
    in_nodes[0] = node(0.0, 0.0, 2.6);
    in_nodes[1] = node(3.0, 0.0, -0.05);

    beam_element my_beam(in_nodes);
    my_beam.print_info();
}
