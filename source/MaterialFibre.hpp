/**
 * @file MaterialFibre.hpp
 * @brief the fibre that is created by the \ref BeamColumnFiberSection.
 * 
 */
#ifndef MATERIAL_FIBRE_HPP
#define MATERIAL_FIBRE_HPP

#include "maths_defaults.hpp"
#include "Material1D.hpp"

/**
 * @brief fibre object that contains a pointer to the material along with fibre area and coordinates.
 * 
 */

class MaterialFibre {
    protected:
        real area; /**< the area of the fibre. */
        real y, z; /**<the x, y coordinates of fibre. */

    public:
        std::unique_ptr<Material1D> material_ptr; /**< a publicly-accessible unique pointer to material object that allows direct calls to fibre material bypassing intermediate interfacing.*/

        /**
         * @brief 
         * 
         */
        MaterialFibre() = default;

        MaterialFibre(Material1D* mat, real area, real y, real z = 0)
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
        void initialise_fibre(Material1D* mat, real area, real y, real z)
        {
            material_ptr = std::make_unique<Material1D>(*mat);
            this->area = area;
            this->y = y;
            this->z = z;
        }

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
};

#endif 