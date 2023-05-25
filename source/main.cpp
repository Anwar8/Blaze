#include <array>
#include "node.hpp"
#include "beam_element.hpp"

int main () {
    std::array<node, 2> in_nodes;
    in_nodes[0] = node(0.0, 0.0, 0.0);
    in_nodes[1] = node(1.0, 0.0, 0.0);

    beam_element my_beam(in_nodes);
    my_beam.print_info();
}
