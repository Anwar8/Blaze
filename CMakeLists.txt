cmake_minimum_required(VERSION 3.16)

# Project and language definitions
project(Blaze LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
# set(CMAKE_CXX_FLAGS_RELEASE "-O3")

# options for the build:
option(BUILD_TESTS "Build test programs" OFF)
option(BUILD_STATIC_LIBS "Build intermediate libraries and link them to executable statically" OFF)
option(OMP "Build with OpenMP support" OFF)
option(KOKKOS "Build with Kokkos" OFF)

option(INCLUDE_GMSH "gmsh exists and Blaze will use it" OFF)
option(VERBOSE "Blaze operate in verbose mode for debugging" OFF)
option(VERBOSE_SLN "Blaze operate in verbose solution mode for debugging" OFF)
option(VERBOSE_STIFFNESSES "Blaze will post stiffness matrix to output stream" OFF)
option(VERBOSE_NLB "Blaze will output U and dU whenever they are solved for by the basic solver" OFF)
option(LF_VERBOSE "Blaze will print LF and iteration every time it enters a new one" OFF)
option(ELEMENT_VERBOSE "Blaze will print most values and data structures produced by nonlinear plasatic element" OFF)

if(INCLUDE_GMSH)
    set(INCLUDE_GMSH 1)
else(INCLUDE_GMSH)
    set(INCLUDE_GMSH 0)
endif(INCLUDE_GMSH)

if(VERBOSE)
    set(VERBOSE 1)
else(VERBOSE)
    set(VERBOSE 0)
endif(VERBOSE)

if(VERBOSE_SLN)
    set(VERBOSE_SLN 1)
else(VERBOSE_SLN)
    set(VERBOSE_SLN 0)
endif(VERBOSE_SLN)

if(VERBOSE_STIFFNESSES)
    set(VERBOSE_STIFFNESSES 1)
else(VERBOSE_STIFFNESSES)
    set(VERBOSE_STIFFNESSES 0)
endif(VERBOSE_STIFFNESSES)

if(VERBOSE_NLB)
    set(VERBOSE_NLB 1)
else(VERBOSE_NLB)
    set(VERBOSE_NLB 0)
endif(VERBOSE_NLB)

if(LF_VERBOSE)
    set(LF_VERBOSE 1)
else(LF_VERBOSE)
    set(LF_VERBOSE 0)
endif(LF_VERBOSE)

if(ELEMENT_VERBOSE)
    set(ELEMENT_VERBOSE 1)
else(ELEMENT_VERBOSE)
    set(ELEMENT_VERBOSE 0)
endif(ELEMENT_VERBOSE)

configure_file(blaze_config.hpp.in ${CMAKE_CURRENT_SOURCE_DIR}/source/blaze_config.hpp)

if(KOKKOS)
if(OMP)
option(Kokkos_ENABLE_OPENMP "To build the OpenMP backend targeting CPUs" ON)
else(OMP)
option(Kokkos_ENABLE_SERIAL "To build the Serial backend targeting CPUs" ON)
endif(OMP)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/kokkos)
endif(KOKKOS)

add_subdirectory(source)

## main exercutable
add_executable(Blaze)
target_sources(Blaze PRIVATE "source/main.cpp")

target_include_directories(Blaze PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/source)

if(KOKKOS)
message(STATUS "Building with Kokkos")
target_link_libraries (Blaze PRIVATE BlazeModel Kokkos::kokkos)
else(KOKKOS)
if(OMP)
message(STATUS "Building with OpenMP support")
find_package(OpenMP REQUIRED)
target_link_libraries(Blaze PUBLIC BlazeModel OpenMP::OpenMP_CXX)
else(OMP)
message(STATUS "Building without OpenMP support")
target_link_libraries(Blaze PUBLIC BlazeModel)
endif(OMP)
endif(KOKKOS)

if(BUILD_TESTS)
    find_package(GTest REQUIRED)
    add_executable(TestBlaze)
    target_sources(TestBlaze PRIVATE "source/tests/my_test.cpp")
    target_include_directories(TestBlaze PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/source)
    target_link_libraries(TestBlaze PUBLIC BlazeModel GTest::gtest GTest::gtest_main)
endif(BUILD_TESTS)

install(TARGETS Blaze)
if(BUILD_TESTS)
    install(TARGETS TestBlaze)
endif(BUILD_TESTS)