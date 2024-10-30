/**
 * @file ElasticPlasticMaterial.hpp
 * @brief class definition that implements an Elastic-Plastic material.
 */

#ifndef ELASTIC_PLASTIC_MATERIAL_HPP
#define ELASTIC_PLASTIC_MATERIAL_HPP
#include "maths_defaults.hpp"
#include "Material1D.hpp"
#include <cmath>

/**
 * @brief A class for elastic-plastic material with hardening.
 */
class ElasticPlasticMaterial : public Material1D {
protected:
    real beta = 0.0; /**< this is the fraction of strain increment \f$ \Delta \epsilon \f$ that takes the stress to the \f$\bf{current}\f$ yield stress.*/
    bool loading = true; /**< is the material currently being continuously loaded by the strain increment? if material is unloading or reloading then this is false.*/

public:
    ElasticPlasticMaterial() = default;
    ElasticPlasticMaterial(real E, real f, real H)
    {
        initialise_material(E, f, H);
    }

    // Copy constructor
    ElasticPlasticMaterial(const ElasticPlasticMaterial& other) : Material1D(other), beta(other.beta), loading(other.loading)
    {}

    /**
     * @brief initialise the material and all its properties.
     * 
     * @param E Young's modulus.
     * @param f yield strength.
     * @param H hardening parameter.
     */
    virtual void initialise_material(real E, real f, real H)
    {
        reference_E0 = E;
        reference_fy0 = f;
        this->H = H;
        
        starting_E = reference_E0;
        starting_E_t = reference_E0;
        starting_fy = reference_fy0;
        starting_fy_bar = reference_fy0;
        starting_stress = 0.0;
        starting_strain = 0.0;
        starting_plastic_strain = 0.0;
        starting_elastic = true;
        initialise_current_state();
    }

    /**
     * 
     * @brief increments the total strain and calculates the yield function, flow, and hardening.
     * @details this is based on Bhatti's Isotropic hardening algorithm (2006, Pp. 396 - 397)
     * @param d_eps the increment in total strain.
     */
    virtual void increment_strain(real d_eps)
    {
        
        initialise_current_state(); // Copy the starting_ state variables into the current state variables
        real delta_s = E * d_eps; // calculate elastic stress increment
        check_if_loading(delta_s); // evaluate loading vs unloading/reloading states
        real s = stress + delta_s; // estimate stress if entire step is elastic
        evolve_yield_surface(); // perform the hardening calculation
        strain += d_eps; // increment the strain by the strain increment


        if (elastic) // if the material starts off in its iteration as elastic
        {
            // we check if it will remain elastic in this step, first
            eval_yield_function(s); // updates the 'elastic' variable based on incremented stress
            
            if (elastic) // this means the material did not yield in this load-step
            {
                stress = s; // so everything is just elastic - just update the stress and exit this function.
                return;
            } else { // the material did, in fact, yield during this step.
                E_t = E * H/ (E + H); // update the tangent modulus as per Bhatti's similar-triangles approach
                beta = (fy_bar - abs(stress))/(abs(s) - abs(stress)); // calculate the value of \f$\beta\f$ which is non-zero
            }
        } else {
            // we only enter here if we are plastic to begin with
            if (loading) // if material is plastic and we are loading, then we continue being plastic - not very interesting.
            {
                beta = 0;
            } else { // otherwise, material is unloading and is thus elastic again
                elastic = true;
                E_t = E;
                stress = s;
                return;
            }
        }
        stress = stress + beta*delta_s + ((E * H)/(E + H)) * (1 - beta) * d_eps;
        calc_plastic_flow(d_eps); // increment the plastic strain
    }

    virtual void check_if_loading(real delta_s)
    {
        loading = stress*delta_s >= 0;
    }

    /**
     * @brief evaluates the yield function which sets the \ref elastic to True if the material yielded.
     * @details \f$ \abs{s} < \bar{\sigma_y}\f$.
     * @param s is the stress state that is being evaluated for yield.
     */
    virtual void eval_yield_function(real s)
    {
        elastic = abs(s) < fy_bar;
    }


    /**
     * @brief sets the current state variables to equal the starting starting (previous) state.
     * 
     */
    virtual void initialise_current_state()
    {
        E = starting_E;
        E_t = starting_E_t;
        fy = starting_fy;
        fy_bar = starting_fy_bar;
        stress = starting_stress;
        strain = starting_strain;
        plastic_strain = starting_plastic_strain;
        elastic = starting_elastic;
    }

    /**
     * @brief calculates plastic strain rate.
     * 
     */
    virtual void calc_plastic_flow(real d_eps)
    {
        plastic_strain = plastic_strain + ((1 - beta)/(1 + (H/E)))*abs(d_eps);
    }
    
    /**
     * @brief calculates the evolution of the yield surface based on the hardening rule and accumulated plastic strain.
     * 
     */
    virtual void evolve_yield_surface()
    {
        fy_bar = fy + H*plastic_strain;
    }

    /**
     * @brief updates the starting state of the material by making the current state the next starting state.
     * 
     */
    virtual void update_starting_state()
    {
        starting_E = E;
        starting_E_t = E_t;
        starting_fy = fy;
        starting_fy_bar = fy_bar;
        starting_stress = stress;
        starting_strain = strain;
        starting_plastic_strain = plastic_strain;
        starting_elastic = elastic;
    }




    /**
     * @brief evolves the temperature of the material, and re-evaluates the material state.
     * 
     * @param dT 
     */
    virtual void increment_temperature(real dT) 
    {
        return;
    }

    /**
     * @brief Get the current Young's modulus of the material.
     * 
     * @return real The current Young's modulus.
     */
    virtual real get_E() const 
    {
        return E;
    }

    /**
     * @brief Get the current tangent Young's modulus of the material.
     * 
     * @return real The current tangent Young's modulus.
     */
    virtual real get_E_t() const
    {
        return E_t;
    }

    /**
     * @brief Get the current yield strength of the material.
     * 
     * @return real The current yield strength.
     */
    virtual real get_fy() const
    {
        return fy;
    }

    /**
     * @brief Get the current yield stress of the material after hardening/softening.
     * 
     * @return real The current yield stress.
     */
    virtual real get_fy_bar() const
    {
        return fy_bar;
    }

    /**
     * @brief Get the current stress in the material.
     * 
     * @return real The current stress.
     */
    virtual real get_stress() const
    {
        return stress;
    }

    /**
     * @brief Get the current total strain in the material.
     * 
     * @return real The current total strain.
     */
    virtual real get_strain() const
    {
        return strain;
    }

    /**
     * @brief 
     * 
     * @return real The current accumulated plastic strain.
     */
    virtual real get_plastic_strain() const
    {
        return plastic_strain;
    }

    /**
     * @brief 
     * 
     * @return bool True if the material is elastic, false otherwise.
     */
    virtual bool is_elastic() const
    {
        return elastic;
    }

    /**
     * @brief Get the starting Young's modulus of the material.
     * 
     * @return real The starting Young's modulus.
     */
    virtual real get_starting_E() const
    {
        return starting_E;
    }

    /**
     * @brief Get the starting tangent Young's modulus of the material.
     * 
     * @return real The starting tangent Young's modulus.
     */
    virtual real get_starting_E_t() const
    {
        return starting_E_t;
    }

    /**
     * @brief Get the starting yield strength of the material.
     * 
     * @return real The starting yield strength.
     */
    virtual real get_starting_fy() const
    {
        return starting_fy;
    }

    /**
     * @brief Get the starting yield stress of the material after hardening/softening.
     * 
     * @return real The starting yield stress.
     */
    virtual real get_starting_fy_bar() const
    {
        return starting_fy_bar;
    }

    /**
     * @brief Get the starting stress in the material.
     * 
     * @return real The starting stress.
     */
    virtual real get_starting_stress() const
    {
        return starting_stress;
    }

    /**
     * @brief Get the starting total strain in the material.
     * 
     * @return real The starting total strain.
     */
    virtual real get_starting_strain() const
    {
        return starting_strain;
    }

    /**
     * @brief Get the starting accumulated plastic strain in the material.
     * 
     * @return real The starting accumulated plastic strain.
     */
    virtual real get_starting_plastic_strain() const
    {
        return starting_plastic_strain;
    }

    /**
     * @brief Check if the material was elastic at the start of the load step.
     * 
     * @return bool True if the material was elastic, false otherwise.
     */
    virtual bool is_starting_elastic() const
    {
        return starting_elastic;
    }
};

#endif