/**
 * @file LocalAxes.hpp
 * @brief local axes object and its functions.
 * 
 * @copyright Copyright (c) 2024
 * 
 */


#ifndef LOCAL_AXES_HPP
#define LOCAL_AXES_HPP
#include "global_coords.hpp"

class LocalAxes : protected GlobalCoords
{
    protected:
        real angle_w_global;
    public:
    // NonlinearTransform(const BasicOrientation& orient) : base_configuration(orient) {}
       LocalAxes(const coords centroid, GlobalCoords& GlobalCoords)
       {
        this->centroid = centroid;
        
       }
        // void set_centroid(coords centroid) : centroid(centroid) {}
}


#endif