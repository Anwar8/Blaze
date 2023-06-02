#include <iostream>
#include "gmsh.h"
#include "mesh.hpp"

int main(int argc, char** argv) {
    // build the geometry
    gmsh::initialize();
    gmsh::model::add("test_model");
    gmsh::model::geo::addPoint(0,0,0,1.0, 1);
    gmsh::model::geo::addPoint(5.0,0,0,1.0, 2);


    gmsh::model::geo::addLine(1, 2, 1);

    
    gmsh::model::geo::synchronize();

    gmsh::model::addPhysicalGroup(1, {1}, 1, "beam_element");
    
    
    gmsh::model::mesh::generate(1);
    gmsh::write("test.msh");


    return 0;
}
