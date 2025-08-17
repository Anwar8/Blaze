<figure style="text-align:center;">
  <img src="notes/images/Blaze logo_blue.png" alt="Blaze Logo" style="width:50%">
</figure>

## Summary
`Blaze` is a finite element method (FEM) program designed to run in parallel on high performance computing infrastructure. It is intended for large problems in the order of 1e5 to 1e7 degrees of freedom. The codebase is written in C++ and provides an interface for user customisation. 

## Table of Contents
- [Summary](#summary)
- [Table of Contents](#table-of-contents)
- [Requirements and dependencies](#requirements-and-dependencies)
- [Instructions for building dependencies on `Cirrus`](#instructions-for-building-dependencies-on-cirrus)
  - [Eigen](#eigen)
  - [Trilinos](#trilinos)
  - [Google Test](#google-test)
- [Building `Blaze`](#building-blaze)
  - [`CMake` flags](#cmake-flags)
  - [Predefined configurations](#predefined-configurations)
  - [Building the documentation](#building-the-documentation)
- [Running `Blaze`](#running-blaze)
- [Testing `Blaze`](#testing-blaze)
- [Example build and run](#example-build-and-run)
- [Directories](#directories)
- [Known issues](#known-issues)


## Requirements and dependencies
- Serial `Blaze` requires `Eigen` and `CMake`.
- Building the tests requires `GoogleTest`. 
- Documentation requires `Doxygen` and `mathjax`.
- `Blaze` with shared memory parallelisation requires `Kokkos` built with any *CPU* backend. `Blaze` does not currently work with GPU backends.
- `Blaze` with distributed parallelisation requires `MPI`, and the `Trilinos` packages `Kokkos`, `Teuchos`, `Tpetra`, and `Amesos2`. The `Trilinos` packages require `MKL`. 

**Warning:** All dependecies need to be built with the same compilers used to build `Blaze`. `Kokkos` is not compatible with `ICC` on `Cirrus`, and so the `GCC` available from within the `intel-20.4` modules will be used. The `intel-20.4` modules are needed as they include `cmkl`.

## Instructions for building dependencies on `Cirrus`
### Eigen
Download [`Eigen 3.4.0`](https://gitlab.com/libeigen/eigen/-/releases/3.4.0) and unpack it in the workspace on Cirrus. Please see [official getting started page](https://libeigen.gitlab.io/eigen/docs-nightly/GettingStarted.html) for detailed instructions.
### Trilinos
1. Clone the entire `Trilinos` source code:
```bash
git clone https://github.com/trilinos/Trilinos.git
```
2. Load the following modules:
```bash
module load intel-20.4/compilers
module load intel-20.4/mpi
module load intel-20.4/cmkl
module load cmake
```
3. Create subdirectories within the `Trilinos` directory for building and installation, and move to the build directory:
```bash
mkdir trilinos-build
mkdir trilinos-install
cd trilinos-build
```
4. Build `Kokkos`, `Teuchos`, `Tpetra`, and `Amesos2` with the following commands:
```bash
cmake -DTPL_ENABLE_MPI=ON -DTPL_ENABLE_MKL=ON -DTPL_MKL_INCLUDE_DIRS="${MKLROOT}/include" -DTPL_MKL_LIBRARIES="${MKLROOT}/lib/intel64/libmkl_rt.so" -DTPL_ENABLE_BLAS=ON -DTPL_ENABLE_LAPACK=ON -DTPL_BLAS_LIBRARIES="${MKLROOT}/lib/intel64/libmkl_rt.so" -DTPL_LAPACK_LIBRARIES="${MKLROOT}/lib/intel64/libmkl_rt.so" -DTrilinos_ENABLE_Tpetra=ON -DTrilinos_ENABLE_Amesos2=ON -DCMAKE_INSTALL_PREFIX=../trilinos-install ..
```
**Note**: the build command above is long because it is important to tell `Trilinos` where to find `BLAS` and `LAPACK`, and that they are provided by `MKL`. This command should work on `Cirrus` and build correctly.

5. Install using multiple cores
```bash
make -j8 install
```
**Warning:** Installation with `make install` will work, but with 1 core the installation process takes about 3 hours on the `Cirrus` logon nodes.

6. Export the installation path of the installed `Trilinos` to enable `CMake` to find it when building `Blaze`:
```bash
export Trilinos_DIR="/work/mdisspt/mdisspt/z#######/Trilinos/trilinos-install"
```
Where `"/work/mdisspt/mdisspt/z#######/Trilinos/trilinos-install"` is replaced with the actual path `Trilinos` is installed in.
### Google Test
1. Load the following modules:
```bash
module load intel-20.4/compilers
module load intel-20.4/mpi
module load cmake
```
2. Clone `Googletest` to chosen directory:
```bash
git clone https://github.com/google/googletest.git
```
3. Build and install with `CMake`:
```bash
cmake .. -DBUILD_GMOCK=ON -DCMAKE_INSTALL_PREFIX=/work/mdisspt/mdisspt/z#######/diss/googletest
make
make install
```
Where `/work/mdisspt/mdisspt/z#######/diss/googletest` is replaced with the desired installation directory.

4. Update the primary `CMakeLists.txt` for `Blaze` to update the installation directory of `Google Test` in the command `list(APPEND CMAKE_PREFIX_PATH "/work/mdisspt/mdisspt/z#######/diss/googletest")` just before `find_package(GTest REQUIRED)`.


## Building `Blaze`
### `CMake` flags
| Flag Name             | Description                                                                                   | Values      | Default |
|-----------------------|----------------------------------------------------------------------------------------------|-------------|---------|
| `BUILD_TESTS`           | Build test programs                                                                          | ON / OFF    | OFF     |
| `BUILD_STATIC_LIBS`     | Build intermediate libraries and link them statically                                         | ON / OFF    | OFF     |
| `WITH_MPI`              | Build the MPI-distributed version of Blaze                                                   | ON / OFF    | OFF     |
| `KOKKOS`                | Build with Kokkos - shared memory parallelism; needs to be built with either `OMP` or `THREADS` for non-serial backend                                                                           | ON / OFF    | OFF     |
| `OMP`                   | Build with OpenMP support                                                                    | ON / OFF    | OFF     |
| `THREADS`               | Build the C++ Threads backend for Kokkos                                                               | ON / OFF    | OFF     |
| `INCLUDE_GMSH`          | Build with support for gmsh  it                                                            | ON / OFF    | OFF     |
| `VERBOSE`               | Blaze will operate in verbose mode for debugging                                                 | ON / OFF    | OFF     |
| `VERBOSE_SLN`           | Blaze will log every step in the solution procedure to the output stream                                         | ON / OFF    | OFF     |
| `LF_VERBOSE`            | Blaze will print LF and iteration every time it enters a new one - convenient for difficult problems where convergence is an issue                            | ON / OFF    | OFF     |
| `VERBOSE_NLB`           | Blaze will output U and dU whenever solved for by the basic solver - **Warning:** do NOT turn on for large problems                      | ON / OFF    | OFF     |
| `VERBOSE_STIFFNESSES`   | Blaze will post stiffness matrix to output stream - **Warning:** do NOT turn on for large problems                                           | ON / OFF    | OFF     |
| `ELEMENT_VERBOSE`       | Blaze will print most values and data structures produced by nonlinear plastic element - **Warning:** do NOT turn on for large problems      | ON / OFF    | OFF     |

### Predefined configurations
The shell script `build.sh` is provided to build `Blaze` with some predefined configurations.
- Standard release version without any flags:
```bash
bash build.sh
``` 

- Release version with any additional flags from table above.
```bash
bash build.sh custom --D####=ON
```
Where `--D####` is any of the configuration flags above. For example:
```bash
bash build.sh custom --DWITH_MPI=ON --DLF_VERBOSE=ON
```
Builds the distributed memory version of `Blaze` and turns on `LF_VERBOSE`.

- Release version with tests and any additional flags from table above.
```bash
bash build.sh tests --D####=ON
```
Where `--D####` is any of the configuration flags above. 

- Debug version with `LF_VERBOSE` and `VERBOSE_SLN` turned on, and any additional configuration flags from above:
```bash
bash build.sh debug --D####=ON
```
Where `--D####` is any of the configuration flags above. 

To debug `Tpetra` matrices, after completing the build of `Blaze` and before running any executables please also run:
```bash
export TPETRA_DEBUG=1
```

- Release version with `Gmsh` support
```bash
bash build.sh mesh
```
### Building the documentation
- To build the documentation:
```bash
doxygen Doxyfile
```
This should place all documentation in `doc`.

## Running `Blaze`
Executables are placed in `bin` by default. 
The `main` `Blaze` executable builds a frame structure and applies a uniform distributed load to its beams. 

To run:
```bash
bin/Blaze
```

To run `Blaze` in parallel (requires having build with the `WITH_MPI` flag):
```
mpirun -n N bin/Blaze
```
where `N` is the number of cores to run `Blaze` with.

`Blaze` main executable is equipped with a series of flags to allow customising the frame, the loading, and the materials. These are:
| Flag Name           | Description                                                        |
|---------------------|--------------------------------------------------------------------|
| `--nbays`             | Number of bays in the frame                                        |
| `--nfloors`           | Number of floors in the frame                                      |
| `--flange_divisions`  | Number of divisions in the flange of the I-section                 |
| `--web_divisions`     | Number of divisions in the web of the I-section                    |
| `--elem_type`         | Element type: LinearElastic, NonlinearElastic, or NonlinearPlastic  |
| `--udl`               | Uniformly distributed load (UDL) applied to beams                  |
| `--beam_divisions`    | Number of divisions per beam                                       |
| `--column_divisions`  | Number of divisions per column                                     |
| `--floor_height`      | Height of each floor                                               |
| `--beam_length`       | Length of each beam                                                |
| `--max_LF`            | Maximum load factor                                                |
| `--nsteps`            | Number of load steps to apply the load over                                               |
| `--tolerance`         | Convergence tolerance                                              |
| `--max_iterations`    | Maximum number of iterations allowable per load step                         |
| `--tf`                | Flange thickness of the I-section                                  |
| `--tw`                | Web thickness of the I-section                                     |
| `--b`                 | Flange width of the I-section                                      |
| `--h`                 | Section height of the I-section                                    |
| `--yield_strength`    | Yield strength of the material                                     |
| `--youngs_modulus`    | Young's modulus of the material                                    |
| `--hardening_ratio`   | Hardening ratio of the material                                    |

## Testing `Blaze`
To run unit tests (cannot run when built with the `WITH_MPI` flag):
```bash
bin/UnitTestBlaze
```

To run verification tests (cannot run when built with the `WITH_MPI` flag):
```bash
bin/VerificationTestsBlaze
```

To run distributed memory tests (requires having built with the `WITH_MPI` flag):
```bash
mpirun -n N bin/TestBlazeMPI
```
Where `N` is any number of processes from 1 to 5. When running parallel tests, expect one pass flag per test per core. So, for a successful run on 3 cores, should see the following output:
```bash
[----------] Global test environment tear-down
[==========] 13 tests from 8 test suites ran. (1194 ms total)
[  PASSED  ] 11 tests.
[  SKIPPED ] 2 tests, listed below:
[  SKIPPED ] DistributedModelSimplySuportedUdlPlastic.CheckMidSpanDeflection
[  SKIPPED ] DistributedModelSimplySuportedUdlElastic.CheckMidSpanDeflection
[==========] 13 tests from 8 test suites ran. (1194 ms total)
[  PASSED  ] 13 tests.
[==========] 13 tests from 8 test suites ran. (1194 ms total)
[  PASSED  ] 11 tests.
[  SKIPPED ] 2 tests, listed below:
[  SKIPPED ] DistributedModelSimplySuportedUdlPlastic.CheckMidSpanDeflection
[  SKIPPED ] DistributedModelSimplySuportedUdlElastic.CheckMidSpanDeflection
```
`DistributedModelSimplySuportedUdlElastic` and `DistributedModelSimplySuportedUdlPlastic` are meant to be skipped on all but one core. For a successful run on `N` cores, expect to see 13 tests passing on one core, and 11 passing with 2 skipped on all other cores.

## Example build and run
The following is an example of building and running a `Blaze` frame model including all expected output.
```bash
[z2259894@cirrus-login2 Blaze]$ module load intel-20.4/compilers
Loading intel-20.4/compilers
  Loading requirement: intel-license gcc/10.2.0 intel-20.4/cc intel-20.4/fc
[z2259894@cirrus-login2 Blaze]$ module load intel-20.4/mpi
[z2259894@cirrus-login2 Blaze]$ module load intel-20.4/cmkl
[z2259894@cirrus-login2 Blaze]$ module load cmake
[z2259894@cirrus-login2 Blaze]$ export Trilinos_DIR="/work/mdisspt/mdisspt/z2259894/Trilinos/trilinos-install"
[z2259894@cirrus-login2 Blaze]$ bash build.sh custom -DWITH_MPI=ON
-- The CXX compiler identification is GNU 10.2.0
-- Detecting CXX compiler ABI info
-- Detecting CXX compiler ABI info - done
-- Check for working CXX compiler: /work/y07/shared/cirrus-software/gcc/10.2.0/bin/c++ - skipped
-- Detecting CXX compile features
-- Detecting CXX compile features - done
-- Enabled Kokkos devices: SERIAL
CMake Warning at /work/mdisspt/mdisspt/z2259894/Trilinos/trilinos-install/lib64/cmake/Kokkos/KokkosConfigCommon.cmake:59 (message):
  The installed Kokkos configuration does not support CXX extensions.
  Forcing -DCMAKE_CXX_EXTENSIONS=Off
Call Stack (most recent call first):
  /work/mdisspt/mdisspt/z2259894/Trilinos/trilinos-install/lib64/cmake/Kokkos/KokkosConfig.cmake:43 (include)
  /work/mdisspt/mdisspt/z2259894/Trilinos/trilinos-install/lib64/cmake/TeuchosCore/TeuchosCoreConfig.cmake:148 (include)
  /work/mdisspt/mdisspt/z2259894/Trilinos/trilinos-install/lib64/cmake/Teuchos/TeuchosConfig.cmake:146 (include)
  /work/mdisspt/mdisspt/z2259894/Trilinos/trilinos-install/lib64/cmake/Amesos2/Amesos2Config.cmake:137 (include)
  /work/mdisspt/mdisspt/z2259894/Trilinos/trilinos-install/lib64/cmake/Trilinos/TrilinosConfig.cmake:92 (include)
  source/core/CMakeLists.txt:4 (find_package)


-- Enabled Kokkos devices: SERIAL
-- Found MPI_CXX: /work/y07/shared/cirrus-software/intel/compilers_and_libraries_2020.4.304/linux/mpi/intel64/lib/libmpicxx.so (found version "3.1") 
-- Found MPI: TRUE (found version "3.1")  
-- Enabled Kokkos devices: SERIAL
-- Enabled Kokkos devices: SERIAL
-- Building without OpenMP support
-- Building with MPI support.
-- Enabled Kokkos devices: SERIAL
-- Enabled Kokkos devices: SERIAL
-- Trilinos include dirs: /work/mdisspt/mdisspt/z2259894/Trilinos/trilinos-install/include
-- Trilinos libraries: Amesos2::all_libs;ShyLU_Node::all_libs;ShyLU_NodeTacho::all_libs;Zoltan2Core::all_libs;Pamgen::all_libs;Galeri::all_libs;Xpetra::all_libs;Thyra::all_libs;ThyraTpetraAdapters::all_libs;ThyraCore::all_libs;TrilinosSS::all_libs;Tpetra::all_libs;TpetraCore::all_libs;TpetraTSQR::all_libs;Zoltan::all_libs;RTOp::all_libs;KokkosKernels::all_libs;Teuchos::all_libs;TeuchosKokkosComm::all_libs;TeuchosKokkosCompat::all_libs;TeuchosRemainder::all_libs;TeuchosNumerics::all_libs;TeuchosComm::all_libs;TeuchosParameterList::all_libs;TeuchosParser::all_libs;TeuchosCore::all_libs
-- Configuring done
-- Generating done
-- Build files have been written to: /work/mdisspt/mdisspt/z2259894/diss/Blaze/build
[ 50%] Building CXX object CMakeFiles/Blaze.dir/source/main.cpp.o
[100%] Linking CXX executable Blaze
[100%] Built target Blaze
Install the project...
-- Install configuration: "Release"
-- Installing: /work/mdisspt/mdisspt/z2259894/diss/Blaze/bin/Blaze
-- Set runtime path of "/work/mdisspt/mdisspt/z2259894/diss/Blaze/bin/Blaze" to ""
[z2259894@cirrus-login2 Blaze]$ mpirun -n 3 bin/Blaze --nbays 4 --nfloors 3 --nsteps 10 
nbays,nfloors,beam_length,floor_height,beam_divisions,column_divisions,element_type,tf,tw,b,h,yield_strength,hardening_ratio,youngs_modulus,flange_divisions,web_divisions,udl,max_LF,nsteps,tolerance,max_iterations,num_nodes,num_elements
4,3,5,3.5,50,35,0,0.0196,0.0114,0.1928,0.4672,4.55e+08,0.01,2e+11,10,40,-3000,1,10,0.01,10,1118,1125

---<Analysis complete. LF = 1, and out-of-balance = 6.73522e-07>---
ANALYSIS_SUCCEEDED
num_iterations:20
rank,mesh_setup,bc_load_records,initialisation,solution,all
0,0.015676975,0.0034379959,0.0087530613,0.09206295,0.12008786
1,0.015559912,0.004144907,0.0083069801,0.091953039,0.11998296
2,0.01555109,0.0033371449,0.0091109276,0.091959,0.11997199
rank,U_to_nodes_mapping,element_state_update,assembly,convergence_check,dU_calculation,material_state_update,result_recording,all
0,0.0014255047,0.020549059,0.03670311,0.00034046173,0.032682419,2.4795532e-05,0.00015926361,0.091956139
1,0.0013539791,0.020038843,0.037235022,0.00042581558,0.032634735,2.3841858e-05,0.000177145,0.091950893
2,0.0013904572,0.02014637,0.037105322,0.00044822693,0.032651424,2.4318695e-05,0.00013685226,0.091958046
```

## Directories
- `./doc` contains the html documentation generated using `Doxygen`. Last updated on 17 August 2025. In future dates, please see [Building the documentation](#building-the-documentation) in [Building `Blaze`](#building-blaze).
- `./bin` contains binary. Automatically cleared when building with `CMake`.
- `./build` intermediate build directory. Auto-cleared by `CMake`.
- `./source` source code and deprecated `Makefile` and `config.mk`.
- `./source/mesh` contains `gmsh` code and `Makefile` for generating a `.msh` file.
- `./source/tests` contains unit tests source code for `Blaze`.

## Known issues
- Number of used threads cannot be output with `read_parallelism_information` when using the `C++ Threads` backend for `Kokkos`. This returns a 0.
- An error where the solver is unable to factorise matrix arises for the Plastic Cantilever test where it takes multiple runs to correctly proceed. This bug is difficult to replicate, has not been explicitly squished, but has not presented recently.



