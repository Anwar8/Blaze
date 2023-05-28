include config.mk
OBJS=$(SRCS:.cpp=.o)
TEST_OBJS=$(TEST_SRCS:.cpp=.o)

## all         : makes the TARGET.
.PHONY : all
all : $(TARGET) $(TEST_TARGET)

## test        : makes the TEST_TARGET in the current directory.
.PHONY : test
test : $(TEST_TARGET)

$(TARGET) : $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS)

$(TEST_TARGET) : $(TEST_OBJS) $(TARGET)
	$(CXX) $(CXXFLAGS) -o $(TEST_TARGET) $(TEST_OBJS)

%.o: %.cpp %.hpp	
	$(CXX) $(CXXFLAGS) $(LDLIBS) $(INCLUDE_FLAGS) -c -o $@ $<

$(TEST_DIR)/%.o: $(TEST_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) $(TEST_LDLIBS) $(TEST_INCLUDE_FLAGS) -c -o $@ $<

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