GMSH_PATH=/opt/homebrew/Cellar/gmsh/4.11.1_1
TARGET=mesh
PRODUCT=test.msh
$(TARGET) :  mesh.cpp mesh.hpp
	g++ -std=c++20 -I$(GMSH_PATH)/include -o $(TARGET) mesh.cpp -L$(GMSH_PATH)/lib -lgmsh

mesh_analysis :  mesh_analysis.cpp
	g++ -std=c++20 -I$(GMSH_PATH)/include -o $(TARGET)_analysis mesh_analysis.cpp -L$(GMSH_PATH)/lib -lgmsh

.PHONY : clean
clean :
	rm -f $(TARGET) $(PRODUCT) $(TARGET)_analysis