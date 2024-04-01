/**
 * @file NonlinearTransform.hpp
 * @brief Nonlinear transformation matrix for geometrically nonlinear analysis.
 * 
 */

#ifndef NONLINEAR_TRANSFORM_HPP
#define NONLINEAR_TRANSFORM_HPP
#include "basic_orientation.hpp"
/**
 * @brief NonlinearTransform object is responsible for creating the geometerically-nonlinear transformation matrix for beam-column analysis.
 * 
 * @details based on Felippa's Nonlinear Finite Element Method
 * @todo should inherit many of the functionality of the Orientation object which only needs to host an axis system and do some simple trig operations.
 */
class NonlinearTransform
{
private:
    BasicOrientation base_configuration;
    BasicOrientation corotated_configuration;
    BasicOrientation current_configuration;
    real Psi; /**< angle \f$\Psi\f$ between base configuration and current (and equally corotated) configuration*/
    real psi; /**< angle \f$\psi\f$ between global axes and base configuration*/
    real phi; /**< summation of \ref Psi and \ref psi - total angle between global axes and current/corotated configuration*/
public:
    /**
     * @brief Construct a new Nonlinear Transform object.
     * 
     * @param orient orientation object for beam-column element.
     */
    NonlinearTransform(const BasicOrientation& orient) : base_configuration(orient) {}
    ~NonlinearTransform();
};


#endif