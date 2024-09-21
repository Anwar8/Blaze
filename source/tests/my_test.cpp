#include "gtest/gtest.h"
#include "ElasticPlasticMaterialTests.hpp"
#include "FibreSectionTests.hpp"
#include "BeamElementTests.hpp"
#include "NodalLoadTests.hpp"
#include "ScribeTest.hpp"
#include "ModelTests.hpp"

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}