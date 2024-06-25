rm -r build
rm -r bin
mkdir build bin bin/mesh
cp  source/mesh/test.msh bin/mesh
cmake -B build -S . -DCMAKE_INSTALL_PREFIX=. -DCMAKE_BUILD_TYPE=Debug
cd build
make install
cd ../bin

