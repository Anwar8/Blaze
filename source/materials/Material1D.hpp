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
    /**
     * @name reference material properties
     * @brief material properties that are initialised and do not change.
     */
    //@{
    real reference_E0; /**< Initial Young's modulus of the material regardless of temperature or plasticity.*/
    real reference_fy0; /**< Initial yield stress of the material.*/
    real H; /**< The \f$ H \f$ hardening parameter for simple kinematic or isotropic hardening.*/ 
    
    //@}
    /**
     * @name state variables
     * @brief material parameters that are updated during state update.
     * @details to avoid having fictitious plastic strain accumulation, we always start with "previous" state - since this
     * is always the state we start from, Blaze calls it the "starting" state. All state members not prefixed are current
     * state variables.
     */
    //@{
    real E;  /**< Current Young's modulus of the material.*/ 
    real E_t; /**< Current tangent Young's modulus of the material.*/ 
    real fy; /**< Current yield strength.*/
    real fy_bar; /**< Current yield stress of the material after hardening/softening - \f$\bar{\sigma}_y\f$ from Bhatti.*/
    real stress; /**< Current stress in the material.*/
    real strain; /**< Current total strain in the material.*/
    real plastic_strain; /**< Current accumulated plastic strain.*/
    bool elastic; /**< Whether the material is currently still elastic or not.*/

    real starting_E; /**< Young's modulus of the material at the start of the load step.*/ 
    real starting_E_t; /**< Tangent Young's modulus of the material at the start of the load step.*/ 
    real starting_fy; /**< yield strength at the start of the load step.*/
    real starting_fy_bar; /**< Yield stress of the material after hardening/softening at the start of the load step.*/
    real starting_stress; /**< Stress in the material at the start of the load step.*/
    real starting_strain; /**< Total strain in the material at the start of the load step.*/
    real starting_plastic_strain; /**< Accumulated plastic strain at the start of the load step.*/
    bool starting_elastic; /**< Whether the material is still elastic or not at the start of the load step.*/
    //@}
    /**
     * @name thermo-mechanical variables
     * @brief currently no implemented.
     */
    real temperature; /**< Current temperature of the material. In degrees Celsius. Default is \f$ 20 ^\circ C\f$.*/
    real thermal_strain; /**< the strain due to thermal expansion \f$\varepsilon_{th} = \alpha \Delta T\f$.*/
    real alpha; /**< Thermal expansion coefficient.*/

public:
    
    virtual ~Material1D() = default;
    Material1D() = default;
    /**
     * @brief initialise the material and all its properties.
     * 
     * @param E Young's modulus.
     * @param f yield strength.
     * @param H hardening parameter.
     */
    virtual void initialise_material(real E, real f, real H) = 0;

    /**
     * @brief increments the total strain and calculates the yield function, flow, and hardening.
     * 
     * @param d_eps the increment in total strain.
     */
    virtual void increment_strain(real d_eps) = 0;

    /**
     * @brief evaluates the yield function which sets the \ref elastic to True if \f$f(\sigma) < 0\f$ and false otherwise.
     * @param s is the stress state that is being evaluated for yield.
     */
    virtual void eval_yield_function(real s) = 0;

    /**
     * @brief calculates plastic strain rate.
     * @param d_eps the increment in strain.
     */
    virtual void calc_plastic_flow(real d_eps) = 0;
    
    /**
     * @brief calculates the evolution of the yield surface based on the hardening rule and accumulated plastic strain.
     * 
     */
    virtual void evolve_yield_surface() = 0;

    /**
     * @brief updates the starting state of the material by making the current state the next starting state.
     * 
     */
    virtual void update_starting_state() = 0;


    /**
     * @brief evolves the temperature of the material, and re-evaluates the material state.
     * 
     * @param dT 
     */
    virtual void increment_temperature(real dT) = 0;

    /**
     * @brief Get the current Young's modulus of the material.
     * 
     * @return real The current Young's modulus.
     */
    virtual real get_E() const = 0;

    /**
     * @brief Get the current tangent Young's modulus of the material.
     * 
     * @return real The current tangent Young's modulus.
     */
    virtual real get_E_t() const = 0;

    /**
     * @brief Get the current yield strength of the material.
     * 
     * @return real The current yield strength.
     */
    virtual real get_fy() const = 0;

    /**
     * @brief Get the current yield stress of the material after hardening/softening.
     * 
     * @return real The current yield stress.
     */
    virtual real get_fy_bar() const = 0;

    /**
     * @brief Get the current stress in the material.
     * 
     * @return real The current stress.
     */
    virtual real get_stress() const = 0;

    /**
     * @brief Get the current total strain in the material.
     * 
     * @return real The current total strain.
     */
    virtual real get_strain() const = 0;

    /**
     * @brief Get the current accumulated plastic strain in the material.
     * 
     * @return real The current accumulated plastic strain.
     */
    virtual real get_plastic_strain() const = 0;

    /**
     * @brief Check if the material is currently elastic or not.
     * 
     * @return bool True if the material is elastic, false otherwise.
     */
    virtual bool is_elastic() const = 0;

    /**
     * @brief Get the starting Young's modulus of the material.
     * 
     * @return real The starting Young's modulus.
     */
    virtual real get_starting_E() const = 0;

    /**
     * @brief Get the starting tangent Young's modulus of the material.
     * 
     * @return real The starting tangent Young's modulus.
     */
    virtual real get_starting_E_t() const = 0;

    /**
     * @brief Get the starting yield strength of the material.
     * 
     * @return real The starting yield strength.
     */
    virtual real get_starting_fy() const = 0;

    /**
     * @brief Get the starting yield stress of the material after hardening/softening.
     * 
     * @return real The starting yield stress.
     */
    virtual real get_starting_fy_bar() const = 0;

    /**
     * @brief Get the starting stress in the material.
     * 
     * @return real The starting stress.
     */
    virtual real get_starting_stress() const = 0;

    /**
     * @brief Get the starting total strain in the material.
     * 
     * @return real The starting total strain.
     */
    virtual real get_starting_strain() const = 0;

    /**
     * @brief Get the starting accumulated plastic strain in the material.
     * 
     * @return real The starting accumulated plastic strain.
     */
    virtual real get_starting_plastic_strain() const = 0;

    /**
     * @brief Check if the material was elastic at the start of the load step.
     * 
     * @return bool True if the material was elastic, false otherwise.
     */
    virtual bool is_starting_elastic() const = 0;

};

#endif // MATERIAL1D_HPP
