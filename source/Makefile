include config.mk
OBJS=$(SRCS:.cpp=.o) 
MSH_OBJS=$(MSH_SRCS:.cpp=.o) 
TEST_OBJS=$(TEST_SRCS:.cpp=.o)

## all         : makes the TARGET and TEST_TARGET.
.PHONY : all
all : $(TARGET) $(TEST_TARGET)

## test        : makes the TEST_TARGET in the current directory.
.PHONY : test
test : $(TEST_TARGET)

$(TARGET) : maths_defaults.o $(OBJS) $(MSH_OBJS)
	@echo "$'\nMaking $(TARGET):"
	@echo ---------------------------------------------
	$(CXX) $(CXXFLAGS) -o $(TARGET) maths_defaults.o $(OBJS) $(MSH_OBJS) $(MESH_LDLIBS)
	@echo ---------------------------------------------

$(TEST_TARGET) : $(TEST_SRCS) $(SRCS) $(filter-out main.cpp, $(MSH_SRCS))
	@echo "$'\nMaking $(TEST_TARGET):"
	@echo ---------------------------------------------
	# rm source/main.o
	$(TEST_CXX) $(CXXFLAGS) $(TEST_INCLUDE_FLAGS) -o $(TEST_TARGET) $(TEST_SRCS) maths_defaults.o $(OBJS) $(filter-out main.o, $(MSH_OBJS)) $(TEST_LDLIBS) $(MESH_LDLIBS)
# $(TEST_CXX) $(CXXFLAGS) $(TEST_INCLUDE_FLAGS) -o $(TEST_TARGET) $(TEST_SRCS) $(SRCS) $(filter-out main.cpp, $(MSH_SRCS)) $(TEST_LDLIBS)
	@echo ---------------------------------------------

%.o: %.cpp %.hpp	
	$(CXX) $(CXXFLAGS) $(INCLUDE_FLAGS) -c -o $@ $<
# $(CXX) $(CXXFLAGS) $(INCLUDE_FLAGS) $(MESH_FLAGS) -c -o $@ $< $(LDLIBS)

maths_defaults.o: maths_defaults.cpp maths_defaults.hpp
	$(CXX) $(CXXFLAGS) $(INCLUDE_FLAGS) -c -o $@ $< 

main.o: main.cpp main.hpp
	$(CXX) $(CXXFLAGS) $(INCLUDE_FLAGS) $(MESH_FLAGS) -c -o $@ $< $(LDLIBS)

global_mesh.o: global_mesh.cpp global_mesh.hpp
	$(CXX) $(CXXFLAGS) $(INCLUDE_FLAGS) $(MESH_FLAGS) -c -o $@ $< $(LDLIBS)

assembler.o: assembler.cpp assembler.hpp
	$(CXX) $(CXXFLAGS) $(INCLUDE_FLAGS) $(MESH_FLAGS) -c -o $@ $<

basic_solver.o: basic_solver.cpp basic_solver.hpp
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