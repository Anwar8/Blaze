/**
 * @file Material1D.hpp
 * @brief the base class for 1D material.
 */

#ifndef MATERIAL1D_HPP
#define MATERIAL1D_HPP
#include "maths_defaults.hpp"

/**
 * @brief A base pure virtual class for 1D material. 
 */
class Material1D {
protected:

    real E0; /**< Initial Young's modulus of the material regardless of temperature or plasticity.*/
    real E_t; /**< Tangent Young's modulus of the material; temperature dependent Young's modulus.*/ 
    real fy0; /**< Initial yield stress of the material.*/
    real fy; /**< Current yield stress of the material after hardening/softening and thermal degradation.*/
    
    real stress; /**< Current stress in the material.*/
    
    real strain; /**< Current total strain in the material.*/
    real plastic_strain; /**< Current plastic strain.*/
    real elastic_strain; /**< The amount of elastic strain in the material.*/

    
    real kinematic_hardening_k; /**< The \f$ \kappa \f$ hardening parameter for simple kinematic hardening.*/ 
    
    bool elastic; /**< Whether the material is still elastic or not.*/

    real temperature; /**< Current temperature of the material. In degrees Celsius. Default is \f$ 20 ^\circ C\f$.*/
    real thermal_strain; /**< the strain due to thermal expansion \f$\varepsilon_{th} = \alpha \Delta T\f$.*/
    real alpha; /**< Thermal expansion coefficient.*/

public:
    /**
     * @brief initialise the material and all its properties.
     * 
     * @param E Young's modulus.
     * @param f yield strength.
     * @param k kinematic hardening parameter.
     */
    virtual void initialise_material(real E, real f, real k) = 0;

    /**
     * @brief increments the total strain and calculates the yield function, flow, and hardening.
     * 
     * @param deps the increment in total strain.
     */
    virtual void increment_strain(real deps) = 0;

    /**
     * @brief evaluates the yield function which sets the \ref elastic to True if \f$f(\sigma) < 0\f$ and false otherwise.
     * 
     */
    virtual void eval_yield_function() = 0;

    /**
     * @brief calculates plastic strain rate.
     * 
     */
    virtual void calc_plastic_flow() = 0;
    
    /**
     * @brief calculates the evolution of the yield surface based on the hardening rule and accumulated plastic strain.
     * 
     */
    virtual void evolve_yield_surface() = 0;


    /**
     * @brief evolves the temperature of the material, and re-evaluates the material state.
     * 
     * @param dT 
     */
    virtual void increment_temperature(real dT) = 0;
};

#endif // MATERIAL1D_HPP
