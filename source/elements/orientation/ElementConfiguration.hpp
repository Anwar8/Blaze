/**
 * @file ElementConfiguration.hpp
 * @brief defines the element configuration class.
 * 
 */

#ifndef ELEMENT_CONFIGURATION_HPP
#define ELEMENT_CONFIGURATION_HPP
#include "maths_defaults.hpp"


/**
 * @brief defines element configuration in terms of end-node coordinates, local axis orientation, and angles with global axis.
 * 
 */
class ElementConfiguration
{
    protected:
    coords centroid = {0, 0, 0}; /**< centroid of the configuration.*/
    coords x_axis, y_axis, z_axis; /**< unit axes for the configuration local axes.*/
    coords pt1, pt2; /**< coordinates of the first and second nodes of the element.*/
    coords pt21; /**< the vector that denotes the difference between the coordinates of \ref pt2 and \ref pt1*/
    real X21, Y21, Z21; /**< the distance between the second and the first points (nodes) for x, y, and z.*/
    real L; /**< the length of the element between its two nodes. Found from \f$ L = \sqrt{X_{21}^2 + Y_{21}^2 + Z_{21}^2}\f$.*/

    real alpha; /**< Angle between the element configuration and the global configuration.*/
    
    public:
    friend class NonlinearTransform;

    /**
     * @brief Just like the constructor, except it is called on an already established object. 
     * 
     * @todo consider whether I should just use the constructor again.
     * 
     * @param node1 
     * @param node2 
     */
    void update_pts(coords node1, coords node2)
    {
        pt1 = node1;
        pt2 = node2;
        calc_distances();
        calc_axes(); 
    }
    /**
     * @brief calculates \ref X21, \ref Y21, \ref Z21, and \ref L from embedded end-point coordinates \ref pt1 and \ref pt2.
     * @details Uses `Eigen3` norm function to calculate L.
     */
    void calc_distances() 
    {
        pt21 = pt2 - pt1;
        X21 = pt21(0);
        Y21 = pt21(1);
        Z21 = pt21(2);

        centroid = pt1 + 0.5*pt21;

        // Uses Eigen3 norm function to retrieve length.
        L = pt21.norm();
        
    }
    /**
     * @brief calculates the unit vectors \ref x_axis, \ref y_axis, and \ref z_axis.
     * @details based on sensible-looking: https://stackoverflow.com/questions/14607640/rotating-a-vector-in-3d-space
     * Rotation around z (to get y-axis from x-axis):
     * 
     * \f$ \boldsymbol{X}^{zz=90 ^{\circ}} = \boldsymbol{T}\boldsymbol{X} = \begin{bmatrix} 0 & -1 & 0 \\ 1 & 0 & 0\\ 0 & 0 & 1\end{bmatrix} \begin{bmatrix} x \\ y \\ z \end{bmatrix} = \begin{bmatrix} -y \\ x \\ z \end{bmatrix}\f$
     *
     * Rotation around y (to get z-axis from x-axis):
     * 
     * \f$ \boldsymbol{X}^{yy=90 ^{\circ}} = \boldsymbol{T}\boldsymbol{X} = \begin{bmatrix} 0 & 0 & 1 \\ 0 & 1 & 0\\ -1 & 0 & 0\end{bmatrix} \begin{bmatrix} x \\ y \\ z \end{bmatrix} = \begin{bmatrix} z \\ y \\ -x \end{bmatrix}\f$
     * 
     * @warning This rotation assumes "Y" is the vertical (gravity) axis. Local axis assumes local z-axis is into the page.
     * @todo make sure this function returns general rotations for 3D. I suspect it already does, but better double-check. Replace with proper reference for transformations.
     */
    void calc_axes()
    {
        x_axis = pt21/L;
        y_axis(0) = -x_axis(1);
        y_axis(1) = x_axis(0);
        y_axis(2) = x_axis(2);

        z_axis(0) = x_axis(2);
        z_axis(1) = x_axis(1);
        z_axis(2) = -x_axis(0);
    }
    

    
};
#endif