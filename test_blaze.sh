cd bin
if [ -z $1 ]; then
    ./TestBlaze
elif [ $1 = "mat" ]; then
	./TestBlaze --gtest_filter="ElasticPlasticMaterialTest.*"
else
    echo "Unknown input: $1. Expected \"mat\" or nothing."
fi