include config.mk
OBJS=$(SRCS:.cpp=.o) 
TEST_OBJS=$(TEST_SRCS:.cpp=.o)

## all         : makes the TARGET.
.PHONY : all
all : $(TARGET) $(TEST_TARGET)

## test        : makes the TEST_TARGET in the current directory.
.PHONY : test
test : $(TEST_TARGET)

$(TARGET) : $(OBJS) global_mesh.o
	@echo "$'\nMaking $(TARGET):"
	@echo ---------------------------------------------
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS) global_mesh.o $(MESH_LDLIBS)
	@echo ---------------------------------------------

$(TEST_TARGET) : $(TEST_SRCS) $(filter-out main.cpp, $(SRCS))
	@echo "$'\nMaking $(TEST_TARGET):"
	@echo ---------------------------------------------
	$(TEST_CXX) $(CXXFLAGS) $(TEST_INCLUDE_FLAGS) -o $(TEST_TARGET) $(TEST_SRCS) $(filter-out main.cpp, $(SRCS)) $(TEST_LDLIBS)
	@echo ---------------------------------------------

%.o: %.cpp %.hpp	
	$(CXX) $(CXXFLAGS) $(INCLUDE_FLAGS) $(MESH_FLAGS) -c -o $@ $< $(LDLIBS)

global_mesh.o: mesh/global_mesh.cpp mesh/global_mesh.hpp
	$(CXX) $(CXXFLAGS) $(INCLUDE_FLAGS) $(MESH_FLAGS) -c -o $@ $<

## clean       : removes all object files and executables for both TARGET and TEST_TARGET.
.PHONY : clean
clean :
	rm -f *.o $(TEST_DIR)/*.o $(TARGET) $(TEST_TARGET)

## clean-test  : removes all object files and executables for TEST_TARGET.
.PHONY : clean-test
clean-test :
	rm -f $(TEST_DIR)/*.o $(TEST_TARGET)

## clean-target: removes all object files and executables for TARGET.
.PHONY : clean-target
clean-target :
	rm -f *.o $(TARGET)

## help        : prints all in-line documentation about current make file.
.PHONY : help
help : Makefile
	@sed -n 's/^##//p' $<