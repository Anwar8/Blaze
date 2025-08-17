rm -r build/CMakeCache.txt

if [ -z $1 ]; then
    cmake -B build -S . -DCMAKE_INSTALL_PREFIX=. -DCMAKE_BUILD_TYPE=Release
elif [ $1 = "custom" ]; then
    shift
    cmake -B build -S . -DCMAKE_INSTALL_PREFIX=. -DCMAKE_BUILD_TYPE=Release "$@"
elif [ $1 = "mesh" ]; then
    cp  source/mesh/test.msh bin/mesh
    cmake -B build -S . -DINCLUDE_GMSH=ON -DCMAKE_INSTALL_PREFIX=. -DCMAKE_BUILD_TYPE=Release
elif [ $1 = "tests" ]; then
    shift
    cmake -B build -S . -DBUILD_TESTS=ON -DCMAKE_INSTALL_PREFIX=. -DCMAKE_BUILD_TYPE=Release "$@"
elif [ $1 = "debug" ]; then    
    shift
    cmake -B build -S . -DBUILD_STATIC_LIBS=OFF -DVERBOSE_SLN=ON -DLF_VERBOSE=ON -DCMAKE_INSTALL_PREFIX=. -DCMAKE_BUILD_TYPE=Debug "$@"
else
    echo "Unknown input: $1. Expected \"custom\", \"mesh\", \"tests\", \"debug", or nothing."
fi
cd build
make install
cd ../bin
