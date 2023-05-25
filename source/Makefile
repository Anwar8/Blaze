CXX = g++
CXXFLAGS = -std=c++2a -I/usr/include/eigen3
SRCS = main.cpp beam_element.cpp node.cpp
OBJS = $(SRCS:.cpp=.o)
TARGET = main

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS)

main.o: main.cpp beam_element.hpp node.hpp
	$(CXX) $(CXXFLAGS) -c main.cpp

beam_element.o: beam_element.cpp beam_element.hpp node.hpp
	$(CXX) $(CXXFLAGS) -c beam_element.cpp

node.o: node.cpp node.hpp
	$(CXX) $(CXXFLAGS) -c node.cpp

clean:
	rm -f $(OBJS) $(TARGET)