#include "gtest/gtest.h"
#include "BeamElementTests.hpp"
#include "NodalLoadTests.hpp"
#include "ScribeTest.hpp"
#include "ModelTests.hpp"
#include "ElasticPlasticMaterialTests.hpp"


int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}