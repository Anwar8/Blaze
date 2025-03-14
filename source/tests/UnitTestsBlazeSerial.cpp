#include "ElasticPlasticMaterialTests.hpp"
#include "FibreSectionTests.hpp"
#include "BeamElementTests.hpp"
#include "NodalLoadTests.hpp"
#include "ScribeTest.hpp"
#include "PlasticBeamElementTests.hpp"
#include "ModelTests.hpp"
#include "PlasticModelTests.hpp"

#include "gtest/gtest.h"

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}