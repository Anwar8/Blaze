cd bin
if [ -z $1 ]; then
    ./TestBlaze
elif [ $1 = "mat" ]; then
	./TestBlaze --gtest_filter="ElasticPlasticMaterialTest.*"
elif [ $1 = "fibre" ]; then
	./TestBlaze --gtest_filter="FibreSectionCentroidTests.*"
    ./TestBlaze --gtest_filter="MaterialFibreTests.*"
else
    echo "Unknown input: $1. Expected \"mat\" or nothing."
fi