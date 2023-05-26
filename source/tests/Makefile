# Makefile for building my_test.cpp with Google Test

# Compiler and linker settings

ifeq ($(shell hostname), Mhds-Air)
CXX = g++
CXXFLAGS = -std=c++20 -I/opt/homebrew/Cellar/eigen/3.4.0_1/include/eigen3 -I$(GTEST_DIR)/include
LDLIBS = -L$(GTEST_DIR)/lib -lgtest -lgtest_main -lpthread
GTEST_DIR = /opt/homebrew/Cellar/googletest/1.13.0
else
CXX = g++
CXXFLAGS = -std=c++2a -I/usr/include/eigen3 -I$(GTEST_DIR)/include
LDLIBS = -L$(GTEST_DIR)/lib -lgtest -lgtest_main -lpthread
GTEST_DIR = /opt/homebrew/Cellar/googletest/1.13.0
endif


# Google Test directory (change as needed)


# Source file and output executable name
SRC = my_test.cpp
EXEC = my_test

# Build rules
all: $(EXEC)

$(EXEC): $(SRC)
	$(CXX) $(CXXFLAGS) $(LDLIBS) $^ -o $@

clean:
	rm -f $(EXEC)