/**
 * @file MaterialFibre.hpp
 * @brief the fibre that is created by the \ref BeamColumnFiberSection.
 * 
 */
#ifndef MATERIAL_FIBRE_HPP
#define MATERIAL_FIBRE_HPP

#include "maths_defaults.hpp"
#include "Material1D.hpp"
#include "ElasticPlasticMaterial.hpp"

/**
 * @brief fibre object that contains a pointer to the material along with fibre area and coordinates.
 * 
 */
class MaterialFibre {
    protected:
        real area; /**< the area of the fibre. */
        real y, z; /**<the x, y coordinates of fibre. */
        real force; /**<the axial force in the fibre. */

    public:
        std::unique_ptr<ElasticPlasticMaterial> material_ptr; /**< a publicly-accessible unique pointer to material object that allows direct calls to fibre material bypassing intermediate interfacing.*/

        /**
         * @brief 
         * 
         */
        MaterialFibre() = default;

        MaterialFibre(ElasticPlasticMaterial* mat, real area, real y, real z = 0)
        {
            initialise_fibre(mat, area, y, z);
        }
        /**
         * @brief used to initialise the fibre object with its properties.
         * 
         * @param mat a pointer to the material object that is copied into this fibre object.
         * @param area the area of the fibre.
         * @param y the y-coordinate of the fibre.
         * @param z the z-coordinate of the fibre.
         */  
        void initialise_fibre(ElasticPlasticMaterial* mat, real area, real y, real z)
        {
            material_ptr = std::make_unique<ElasticPlasticMaterial>(*mat);
            this->area = area;
            this->y = y;
            this->z = z;
        }

        /**
         * @brief a copy-constructor for the fibre that allows us to copy the material object pointed to by the unique_ptr.
         * 
         * @param other The MaterialFibre object to be copied.
         */
        MaterialFibre(const MaterialFibre& other)
        {
            area = other.area;
            y = other.y;
            z = other.z;
            material_ptr = std::make_unique<ElasticPlasticMaterial>(*other.material_ptr);
        }
        
        
        // void increment_material_strain(real strain_increment);
        // void get_material_stress();
        // void get_material_Et();
        // void update_material_starting_state();
        // Material1D* get_ptr_to_material();

        /**
         * @brief Get the area of the fibre.
         * 
         * @return The area of the fibre.
         */
        real get_area() const { return area; }

        /**
         * @brief Get the y-coordinate of the fibre.
         * 
         * @return The y-coordinate of the fibre.
         */
        real get_y() const { return y; }

        /**
         * @brief Get the z-coordinate of the fibre.
         * 
         * @return The z-coordinate of the fibre.
         */
        real get_z() const { return z; }
        
        void calc_force() {
            force = material_ptr->get_stress() * area;
        }
        real get_force() {return force;}

        void print_info()
        {
            std::cout << "Fibre at (y,z) = (" << y << ", " << z << "), with A = " << area << ", and force = " << force << "." << std::endl;
            std::cout << "Its material has E = " << material_ptr->get_E() << ", and fy = " << material_ptr->get_fy() << "." << std::endl;
        }
};

#endif 


