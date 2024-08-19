/**
 * @file BeamColumnFiberSection.hpp
 * @brief definitions for the fibre-based beam-column element sections
 * 
 */
#ifndef BEAM_COLUMN_FIBER_SECTION_HPP
#define BEAM_COLUMN_FIBER_SECTION_HPP

#include "maths_defaults.hpp"
#include "MaterialFibre.hpp"
#include "Material1D.hpp"
/**
 * @brief cross-section containing fibres used to represent material nonlinearity of a beam-column element.
 * 
 */

class BeamColumnFiberSection {
    private:
        real section_area; /**<total area of the section - combined area of all fibres.*/ 
        real weighted_E; /**<an equivalent Young's modulus taken as the weighted mean of all fibres.*/
        
        std::vector<MaterialFibre> fibres; /**<a vector containing all fibres of the cross-section.*/
        
        real moment_yy; /**< the moment of the section about its y axis.*/
        real axial_force; /**< the axial force in the section.*/

        real axial_strain; /**< Current axial strain \f$\varepsilon_{axial}\f$.*/
        real curvature;  /**< Current curvature \f$\kappa\f$.*/
        
        real starting_axial_strain; /**< Current axial strain \f$\varepsilon_{axial}\f$.*/
        real starting_curvature;  /**< Current curvature \f$\kappa\f$.*/

        real y_bar; /**< distance to the centroid of the section relative to the plane at which y = 0.*/

        mat  D_t = make_xd_mat(2,2); /**< the 2x2 tangent constitutive matrix of the section.*/
        
    public:
        /**
         * @brief populates the fibre vector of the section.
         * 
         * @tparam stl_container an STL-compatible container such as a std::vector. Needs to have an interator.
         * @param mat a pointer to a material object that will be copied into each fibre.
         * @param areas the area of each fibre.
         * @param ys the y-coordinate of each fibre in order.
         */
        template <typename stl_container, typename MaterialType>
        void add_fibres(MaterialType* mat, stl_container areas, stl_container ys)
        {
            if (areas.size() != ys.size())
            {
                std::cout << "BeamColumnFiberSection::add_fibres can only take equally-sized input arrays. sizes of area and ys are: " << areas.size() << ", " << ys.size() << "." << std::endl;
                exit(1);
            }
            auto y_iterator = ys.begin();
            
            for (auto area: areas)
            {
                fibres.emplace_back(MaterialFibre(mat, area, *y_iterator));
                ++y_iterator;
            }
        }

        /**
         * @brief increment the section axial strain and curvature.
         * 
         * @param new_axial_strain 
         * @param new_curvature 
         */
        void increment_section_strains(real new_axial_strain, real new_curvature)
        {
            axial_strain = new_axial_strain;
            curvature = new_curvature;
        }

        /**
         * @brief calculates the strain for each fibre and applies it. 
         * @details as per equation (18) from izzuddin's notes (chapter on Materially Nonlinear Finite Elements): 
         * \f$ \boldsymbol{\varepsilon}_{fibre} = \boldsymbol{y}_{fibre} \boldsymbol{\varepsilon} = \begin{bmatrix} \varepsilon_1 \\ \varepsilon_2 \\ \vdots \\ \varepsilon_n\end{bmatrix} =  \begin{bmatrix} 1 & y_1 \\ 1 & y_2 \\ \vdots & \vdots \\ 1 & y_n \end{bmatrix} \begin{bmatrix} \varepsilon_{axial} \\ \kappa \end{bmatrix}\f$
         */
        void increment_fibre_strains()
        {
            real d_axial_strain = axial_strain - starting_axial_strain;
            real d_curvature = curvature - starting_curvature;

            // careful - since each fibre has a unique_ptr, it cannot be copied and so must be accessed by reference even in the loop.
            for (auto& fibre : fibres)
            {
                real strain_increment = d_axial_strain - (fibre.get_y() - y_bar)*d_curvature; // the minus sign is as per Izzuddin's notation.
                fibre.material_ptr->increment_strain(strain_increment);
            }
        }

        /**
         * @brief calculates the axial forces and moment of the section from the contribution of the different fibres. Requires a known section centroid.
         * 
         */
        void calc_section_forces()
        {
            axial_force = 0.0;
            moment_yy = 0.0;
            for (auto& fibre : fibres)
            {
                // std::cout << "fibre y = " << fibre.get_y() << " has area " << fibre.get_area() << " and has stress = " << fibre.material_ptr->get_stress() << std::endl;
                real force = (fibre.material_ptr->get_stress())*fibre.get_area();
                axial_force += force;
                moment_yy += force*-(fibre.get_y() - y_bar);
            }
        }

        /**
         * @brief calculates an area-weighted-mean of the Young's modulus of the section. 
         * @details \f$ E_{section} = \frac{ \sum_{fibre = i}^n A_{i} E_{i}}{A_{section}}\f$
         */
        void calc_area_weighted_E()
        {
            section_area = 0.0;
            real area_times_E = 0.0;
            for (auto& fibre : fibres)
            {
                real fibre_area = fibre.get_area();
                section_area += fibre_area;
                area_times_E += fibre_area*fibre.material_ptr->get_E_t();
            }
            weighted_E = area_times_E/section_area;
        }

        /**
         * @brief calculates the distance to the centroid \ref y_bar of the section.
         * @details This calculation accounts for changes in Young's modulus due to plasticity and temperature: \f$ \bar{y} = \frac{ \sum_{fibre = i}^n A_{i} (E_{t,i}/E_{section})}{A_{section}} = \frac{ \sum_{fibre = i}^n A_{i} E_{t,i}}{A_{section}E_{section}}\f$
         */
        void calc_section_centroid()
        {
            real area_moment = 0.0;
            for (auto& fibre : fibres)
            {
                area_moment += fibre.get_y()*fibre.get_area()*(fibre.material_ptr->get_E_t());
            }
            y_bar = area_moment/(section_area*weighted_E);
        }

        /**
         * @brief calculates the tangent constitutive matrix \ref D_t.
         * @details Equation (20.a) from Izzuddin's notes:
         * 
         * \f$ \boldsymbol{D}_t = \boldsymbol{Y}^T \boldsymbol{A} \boldsymbol{E}_t \boldsymbol{Y} \f$
         * 
         * Which when expanded returns the very simple 2x2 matrix:
         * \f$ \boldsymbol{D}_t = \begin{bmatrix} \sum_i^n A_i E_{t,i} & \sum_i^n -A_i E_{t,i} y_i \\ \sum_i^n -A_i E_{t,i} y_i & \sum_i^n A_i E_{t,i} y_i^2  \end{bmatrix}\f$
         * 
         * With the careful consideration that in Izzuddin's notes, \f$ y_i\f$ is the distance between the fibre and centroid, so we will need to substitute:
         * 
         * \f$ y_i^{(Izzuddin)} = y_i - \bar{y}\f$, so we will rewrite \f$ \boldsymbol{D}_t\f$ as:
         * 
         * \f$ \boldsymbol{D}_t = \begin{bmatrix} \sum_i^n A_i E_{t,i} & \sum_i^n -A_i E_{t,i} (y_i - \bar{y}) \\ \sum_i^n -A_i E_{t,i} (y_i - \bar{y}) & \sum_i^n A_i E_{t,i} (y_i - \bar{y})^2  \end{bmatrix}\f$
         */
        void calc_tan_contitutive_matrix()
        {             
            D_t.setZero();
            for (auto& fibre : fibres)
            {
                real E_t_i = fibre.material_ptr->get_E_t();
                real A_i = fibre.get_area();
                real y_i = fibre.get_y();

                D_t(0,0) = D_t(0,0) + A_i * E_t_i;
                D_t(1,1) = D_t(1,1) + A_i * E_t_i * pow(y_i - y_bar, 2); 
                D_t(1,0) = D_t(1,0) - A_i * E_t_i * (y_i - y_bar); 
                D_t(0,1) = D_t(1,0);
            }
        }

        /**
         * @brief applies updates the section state by applying an a strain vector after calculating the section centroid and incrementing the section strains. 
         * 
         * @param epsilon a 2-row vector containing axial strain and curvature.
         */
        void update_section_state(vec& epsilon)
        {
            calc_area_weighted_E();
            calc_section_centroid();

            increment_section_strains(epsilon(0), epsilon(1));
            increment_fibre_strains();
            calc_section_forces();
            calc_tan_contitutive_matrix();
        }

        /**
         * @brief updates the section starting state as well as the fibre states.
         * @details the strains are calculated at the element level for the displacement state not as an increment. Each section, however,
         * calculates the increment in strains from a starting state of the section, and then applies them to the material of each fibre.
         * This function allows for 'committing' the state - that is, updating the starting strain for which the strain increments are calcualted.
         */
        void update_section_starting_state()
        {
            starting_axial_strain = axial_strain;
            starting_curvature = curvature;

            for (auto& fibre: fibres)
            {
                fibre.material_ptr->update_starting_state();
            }
        }
        /**
         * @brief Get the total area of the section.
         * 
         * @return The total area of the section.
         */
        real get_section_area() const { return section_area; }

        /**
         * @brief Get the equivalent Young's modulus of the section.
         * 
         * @return The equivalent Young's modulus of the section.
         */
        real get_weighted_E() const { return weighted_E; }

        /**
         * @brief Get the moment of the section about its y axis.
         * 
         * @return The moment of the section about its y axis.
         */
        real get_moment_yy() const { return moment_yy; }

        /**
         * @brief Get the axial force in the section.
         * 
         * @return The axial force in the section.
         */
        real get_axial_force() const { return axial_force; }

        /**
         * @brief Get the current axial strain.
         * 
         * @return The current axial strain.
         */
        real get_axial_strain() const { return axial_strain; }

        /**
         * @brief Get the current curvature.
         * 
         * @return The current curvature.
         */
        real get_curvature() const { return curvature; }

        /**
         * @brief Get the starting axial strain.
         * 
         * @return The starting axial strain.
         */
        real get_starting_axial_strain() const { return starting_axial_strain; }

        /**
         * @brief Get the starting curvature.
         * 
         * @return The starting curvature.
         */
        real get_starting_curvature() const { return starting_curvature; }

        /**
         * @brief Get the distance to the centroid of the section relative to the y = 0 plane.
         * 
         * @return The distance to the centroid of the section relative to the y = 0 plane.
         */
        real get_y_bar() const { return y_bar; }

        /**
         * @brief Get the tangent constitutive matrix of the section.
         * 
         * @return The tangent constitutive matrix of the section.
         */
        const mat& get_D_t() const { return D_t; }
};
#endif