rm -r build
rm -r bin
mkdir build bin bin/mesh
cp  source/mesh/test.msh bin/mesh
cmake -B build -S . -DCMAKE_INSTALL_PREFIX=.
cd build
make install

