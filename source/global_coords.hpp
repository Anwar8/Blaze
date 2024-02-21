/**
 * @file global_coords.hpp
 * @brief centroid and the global unit vectors for x, y, and z
 * 
 */
#ifndef GLOBAL_COORDS
#define GLOBAL_COORDS

#include "maths_defaults.hpp"
/**
 * @brief a class to contain the global coordinate system using unit vectors.
 * 
 * @details The global coordinate system can allow members to transform
 * themselves with respect to it using their /ref BasicOrientation object.
 * 
 */
class GlobalCoords {
    private:
        coords centroid = {0.0, 0.0 , 0.0};
        coords unit_x = {1.0, 0.0, 0.0};
        coords unit_y = {0.0, 1.0, 0.0};
        coords unit_z = {0.0, 0.0, 1.0};
    public: 
        coords get_centroid() {return centroid;}
        coords get_unit_x() {return unit_x;}
        coords get_unit_y() {return unit_y;}
        coords get_unit_z() {return unit_z;}
};
#endif 