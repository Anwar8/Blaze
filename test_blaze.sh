cd bin
if [ -z $1 ]; then
    ./TestBlaze
elif [ $1 = "beam" ]; then
	./TestBlaze --gtest_filter="RigidBodyMotionTest.*"
    ./TestBlaze --gtest_filter="BasicTransformationTest.*"
    ./TestBlaze --gtest_filter="ConstantStrainStateTest.*"
    ./TestBlaze --gtest_filter="ElementMappingTest.*"
elif [ $1 = "mat" ]; then
	./TestBlaze --gtest_filter="ElasticPlasticMaterialTest.*"
elif [ $1 = "fibre" ]; then
	./TestBlaze --gtest_filter="FibreSectionCentroidTests.*"
    ./TestBlaze --gtest_filter="MaterialFibreTests.*"
    ./TestBlaze --gtest_filter="FibreSectionPureBendingTests.*"
    ./TestBlaze --gtest_filter="FibreSectionPureAxialTests.*"
    ./TestBlaze --gtest_filter="FibreSectionTangentConstitutiveMatrix.*"
    ./TestBlaze --gtest_filter="FibreSectionTangentConstitutiveMatrixHardening.*"
elif [ $1 = "plastic" ]; then
	./TestBlaze --gtest_filter="PlasticBeamTests.*"
elif [ $1 = "modelp" ]; then
	./TestBlaze --gtest_filter="MeshTestsPlastic.*"
    ./TestBlaze --gtest_filter="CantileverBeamPlastic.*"
    ./TestBlaze --gtest_filter="SimplySupportedPlastic.*"
    ./TestBlaze --gtest_filter="SimplySupportedUdlPlastic.*"
elif [ $1 = "cant" ]; then
    ./TestBlaze --gtest_filter="CantileverBeamPlastic.*"
else
    echo "Unknown input: $1. Expected \"mat\" or nothing."
fi