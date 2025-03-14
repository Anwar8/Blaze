cd bin
if [ -z $1 ]; then
    ./UnitTestBlaze
elif [ $1 = "beam" ]; then
	./UnitTestBlaze --gtest_filter="RigidBodyMotionTest.*"
    ./UnitTestBlaze --gtest_filter="BasicTransformationTest.*"
    ./UnitTestBlaze --gtest_filter="ConstantStrainStateTest.*"
    ./UnitTestBlaze --gtest_filter="ElementMappingTest.*"
elif [ $1 = "mat" ]; then
	./UnitTestBlaze --gtest_filter="ElasticPlasticMaterialTest.*"
elif [ $1 = "fibre" ]; then
	./UnitTestBlaze --gtest_filter="FibreSectionCentroidTests.*"
    ./UnitTestBlaze --gtest_filter="MaterialFibreTests.*"
    ./UnitTestBlaze --gtest_filter="FibreSectionPureBendingTests.*"
    ./UnitTestBlaze --gtest_filter="FibreSectionPureAxialTests.*"
    ./UnitTestBlaze --gtest_filter="FibreSectionTangentConstitutiveMatrix.*"
    ./UnitTestBlaze --gtest_filter="FibreSectionTangentConstitutiveMatrixHardening.*"
elif [ $1 = "plastic" ]; then
	./UnitTestBlaze --gtest_filter="PlasticBeamTests.*"
elif [ $1 = "modelp" ]; then
	./UnitTestBlaze --gtest_filter="MeshTestsPlastic.*"
    ./VerificationTestsBlaze --gtest_filter="CantileverBeamPlastic.*"
    ./VerificationTestsBlaze --gtest_filter="SimplySupportedPlastic.*"
    ./VerificationTestsBlaze --gtest_filter="SimplySupportedUdlPlastic.*"
elif [ $1 = "cant" ]; then
    ./VerificationTestsBlaze --gtest_filter="CantileverBeamPlastic.*"
elif [ $1 = "unit_custom" ]; then
    shift
    ./UnitTestBlaze --gtest_filter="$@"
elif [ $1 = "verification_custom" ]; then
    shift
    ./VerificationTestsBlaze --gtest_filter="$@"
elif [ $1 = "verification_custom" ]; then
    shift
    ./DistributedTestsBlaze --gtest_filter="$@"
elif [ $1 = "mpi" ]; then
    mpirun -n $2 ./TestBlazeMPI
else
    echo "Unknown input: $1. Expected \"mat\", \"beam\", \"fibre\", \"fibre\", \"plastic\", \"modelp\", \"cant\", \"custom\", or nothing."
fi