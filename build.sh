rm -r build
rm -r bin

if [ -z $1 ]; then
    mkdir build bin
    cmake -B build -S . -DCMAKE_INSTALL_PREFIX=. -DCMAKE_BUILD_TYPE=Release
elif [ $1 = "mesh" ]; then
    mkdir build bin bin/mesh
    cp  source/mesh/test.msh bin/mesh
    cmake -B build -S . -DINCLUDE_GMSH=ON -DCMAKE_INSTALL_PREFIX=. -DCMAKE_BUILD_TYPE=Release 
elif [ $1 = "tests" ]; then
    mkdir build bin
    cmake -B build -S . -DBUILD_TESTS=ON -DCMAKE_INSTALL_PREFIX=. -DCMAKE_BUILD_TYPE=Release 
elif [ $1 = "all" ]; then
    mkdir build bin bin/mesh
    cp  source/mesh/test.msh bin/mesh
    cmake -B build -S . -DINCLUDE_GMSH=ON  -DBUILD_TESTS=ON -DCMAKE_INSTALL_PREFIX=. -DCMAKE_BUILD_TYPE=Release 
else
    echo "Unknown input: $1. Expected \"mesh\", \"tests\" \"all\", or nothing."
fi
cd build
make install
cd ../bin