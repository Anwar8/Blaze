cmake_minimum_required(VERSION 3.16)

# Project and language definitions
project(Blaze LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# options for the build:
option(INCLUDE_GMSH "gmsh exists and Blaze will use it" OFF)
option(BUILD_TESTS "Build test programs" OFF)

# # Getting external packages
# ## Eigen for numerical methods
# find_package(Eigen3 REQUIRED)
# ## gmsh library and header files
# if(INCLUDE_GMSH)
#     find_library(GMSH_LIB Gmsh PATH_SUFFIXES lib)
#     find_path(GMSH_INC "Gmsh.h" PATH_SUFFIXES include */include gmsh include/gmsh)
# endif(INCLUDE_GMSH)

add_subdirectory(source)


## defining the internal library BLAZE_LIB which contains the object files for all source files
# --------------------------------------------------------------------------------------------------------------
### BLAZE_LIB will need to link to the libraries in the subdirectories: BlazeCore BlazeMaterials BlazeElements BlazeUtils; in fact, everything here will need to link to BlazeUtils for the maths
# --------------------------------------------------------------------------------------------------------------
# add_library(BLAZE_LIB STATIC ${SOURCES})
# if(INCLUDE_GMSH)
#     target_include_directories(BLAZE_LIB PUBLIC "source" ${GMSH_INC})
#     target_link_libraries(BLAZE_LIB Eigen3::Eigen ${GMSH_LIB})
# else(INCLUDE_GMSH)
#     target_include_directories(BLAZE_LIB PUBLIC "source")
#     target_link_libraries(BLAZE_LIB Eigen3::Eigen)
# endif(INCLUDE_GMSH)

# Defining the project executables 
## main exercutable
add_executable(Blaze)
target_sources(Blaze PRIVATE "source/main.cpp")

target_link_libraries(Blaze PUBLIC BlazeModel)
target_include_directories(Blaze PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/source)
# if(INCLUDE_GMSH)
#     target_include_directories(Blaze PUBLIC "source" ${GMSH_INC})
#     target_link_libraries(Blaze BLAZE_LIB Eigen3::Eigen ${GMSH_LIB})
# else(INCLUDE_GMSH)
#     target_include_directories(Blaze PUBLIC "source")
#     target_link_libraries(Blaze BLAZE_LIB Eigen3::Eigen)
# endif(INCLUDE_GMSH)
## testing exectuable

# if(BUILD_TESTS)
#     enable_testing()
#     find_package(GTest REQUIRED)
#     add_executable(TestBlaze)
#     target_sources( TestBlaze PRIVATE "source/tests/my_test.cpp")
#     if(INCLUDE_GMSH)
#         target_include_directories(TestBlaze PUBLIC "source" ${GMSH_INC})
#         target_link_libraries(TestBlaze BLAZE_LIB GTest::gtest GTest::gtest_main Eigen3::Eigen ${GMSH_LIB})
#     else(INCLUDE_GMSH)
#         target_include_directories(TestBlaze PUBLIC "source")
#         target_link_libraries(TestBlaze BLAZE_LIB GTest::gtest GTest::gtest_main Eigen3::Eigen)
#     endif(INCLUDE_GMSH)
#     add_test(AllTests TestBlaze) 
# endif(BUILD_TESTS)


# Install directives - what to install?
# if(BUILD_TESTS)
#     install(TARGETS Blaze TestBlaze)
# else(BUILD_TEST)
#     install(TARGETS Blaze)
# endif(BUILD_TESTS)
install(TARGETS Blaze)