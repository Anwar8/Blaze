# Blaze

## Summary
`Blaze` is a finite element method (FEM) program developed specifically for structures in fire. It is explicitly designed from the ground up for scalability on high performance computing (HPC) facilities. The codebase is written in C++ and provides an interface for user customisation.

## Requirements
`Blaze` requires `gmsh`, `Eigen`, and `GoogleTest`. Building requires `CMake`.
Documentation requires `Doxygen` and `mathjax`.

## Usage
To build:
```
bash build.sh
```
To update binary due to minor code change:
```
bash update.sh
```
To run:
```
cd bin
./Blaze
```
To run tests:
```
cd bin
./Test_Blaze
```


## Directories
- `./doc` contains the html documentation generated using `Doxygen`.
- `./notes` additional `Markdown` files used for project management and record keeping.
- `./bin` contains binary. Automatically cleared when building with `CMake`.
- `./build` intermediate build directory. Auto-cleared by `CMake`.
- `./source` source code and deprecated `Makefile` and `config.mk`.
- `./source/mesh` contains `gmsh` code and `Makefile` for generating a `.msh` file.
- `./source/tests` contains unit tests source code for `Blaze`.






