# Blaze

<figure style="text-align:center;">
  <img src="notes/images/Blaze logo_blue.png" alt="Blaze Logo" style="width:50%">
</figure>

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
To run specific tests:
```
bash test_blaze.sh KEYWORD
``` 
where `KEYWORD` is replaced by:
- nothing, to run all tests.
- "mat" to run `ElasticPlasticMaterialTest`.
- "fibre" to run `FibreSectionCentroidTests`, `MaterialFibreTests`, `FibreSectionPureBendingTests`, and `FibreSectionPureAxialTests`.

You can also run `source aliases.sh` to add `build`, `update` and `test` aliases.

## Directories
- `./doc` contains the html documentation generated using `Doxygen`. Not tracked, so will need to be generated fresh.
- `./notes` additional `Markdown` files used for project management and record keeping.
- `./bin` contains binary. Automatically cleared when building with `CMake`.
- `./build` intermediate build directory. Auto-cleared by `CMake`.
- `./POC` proof of concept folder containing python scripts for checking algorithms and implementation before implementation in `Blaze`.
- `./source` source code and deprecated `Makefile` and `config.mk`.
- `./source/mesh` contains `gmsh` code and `Makefile` for generating a `.msh` file.
- `./source/tests` contains unit tests source code for `Blaze`.






