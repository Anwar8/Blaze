#include <iostream>
#include "gtest/gtest.h"
#include "../beam_element.hpp"
#include "../maths_defaults.hpp"
#include "../node.hpp"


TEST(RigidBodyMotion, TestFunction2) {
    
    std::array<node, 2> in_nodes;
    in_nodes[0] = node(0.0, 0.0, 0.0);
    in_nodes[1] = node(3.0, 0.0, 0.0);

    beam_element my_beam(in_nodes);
    my_beam.calc_N(1.5);
    my_beam.calc_B(1.5);
    vec d = make_xd_vec(6);
    d.setZero();

    d(0) = 1; 
    d(3) = 1;
    my_beam.set_d(d);
    my_beam.calc_eps();
    vec eps = my_beam.get_eps();
    EXPECT_FLOAT_EQ(eps(0), 0.0);
    EXPECT_FLOAT_EQ(eps(1), 0.0);
}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}