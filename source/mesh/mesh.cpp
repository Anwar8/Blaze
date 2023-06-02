#include <iostream>
#include "gmsh.h"
#include "mesh.hpp"

int main(int argc, char** argv) {
    // build the geometry
    gmsh::initialize();
    gmsh::model::add("test_model");
    gmsh::model::geo::addPoint(0,0,0,1.0, 1);
    gmsh::model::geo::addPoint(5.0,0,0,1.0, 2);

    gmsh::model::geo::addPoint(6.0,0,0,1.0, 3);
    gmsh::model::geo::addPoint(5.5,0,1.0,1.0, 4);
    
    
    gmsh::model::geo::addPoint(6.5,0,0.0,1.0, 5);
    gmsh::model::geo::addPoint(6.5,0,-0.5,1.0, 6);
    gmsh::model::geo::addPoint(6.0,0,-0.5,1.0, 7);

    gmsh::model::geo::addLine(1, 2, 1);

    gmsh::model::geo::addLine(2, 3, 2);
    gmsh::model::geo::addLine(3, 4, 3);
    gmsh::model::geo::addLine(4, 2, 4);

    gmsh::model::geo::addLine(3, 5, 5);
    gmsh::model::geo::addLine(5, 6, 6);
    gmsh::model::geo::addLine(6, 7, 7);
    gmsh::model::geo::addLine(7, 3, 8);
    



    
    gmsh::model::geo::addCurveLoop({2,3,4}, 1);
    gmsh::model::geo::addCurveLoop({5,6,7,8}, 2);
    
    gmsh::model::geo::addPlaneSurface({1}, 1);
    gmsh::model::geo::addPlaneSurface({2}, 2);
    
    gmsh::model::geo::synchronize();

    gmsh::model::addPhysicalGroup(1, {1}, 1, "beam_element");
    
    gmsh::model::addPhysicalGroup(2, {1}, 2, "plane_surface_triangle");

    gmsh::model::addPhysicalGroup(2, {2}, 3, "plane_surface_rectangle");
    
    gmsh::model::mesh::generate(2);
    gmsh::write("test.msh");


    return 0;
}
