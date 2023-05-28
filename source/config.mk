TARGET= xblaze
SRCS= maths_defaults.cpp node.cpp beam_element.cpp main.cpp

TEST_TARGET=test_xblaze
TEST_DIR=tests
TEST_SRCS= my_test.cpp


# the := operator is necessary here to ensure the 
# directory is appended to the test files
TEST_SRCS := $(addprefix $(TEST_DIR)/,$(TEST_SRCS))

# Choose paths for libraries and frameworks based on where the code is being compiled
ifeq ($(shell hostname), Mhds-Air)
# include paths
EIGEN_PATH= /opt/homebrew/Cellar/eigen/3.4.0_1/include/eigen3
GTEST_DIR = /opt/homebrew/Cellar/googletest/1.13.0
# compiler and settings
CXX= g++-12
TEST_CXX= g++
CXXFLAGS= -std=c++20
else
# include paths
EIGEN_PATH= /opt/homebrew/Cellar/eigen/3.4.0_1/include/eigen3
GTEST_DIR= 
# compiler and settings
CXX= g++
TEST_CXX= g++
CXXFLAGS= -std=c++2a
endif

INCLUDE_FLAGS= -I$(EIGEN_PATH)
# no linker libraries are needed
LDLIBS= 
TEST_INCLUDE_FLAGS= -I$(GTEST_DIR)/include -I$(EIGEN_PATH) -pthread
TEST_LDLIBS= -lgtest -lgtest_main -L$(GTEST_DIR)/lib 