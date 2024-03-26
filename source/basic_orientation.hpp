/**
 * @file basic_orientation.hpp
 * @brief the basic orientation object and functions for element transformations
 * 
 */
#ifndef BASIC_ORIENTATION
#define BASIC_ORIENTATION
#include "maths_defaults.hpp"
#include "node.hpp"

/**
 * @brief 
 * 
 * creates transform matrices to shift and rotate elements
 * 
 * @details
 * 
 * calculates the local axis of a beam-column element, and then creates a transform
 * matrix T that will be used to transform the stiffness matrix.
 * 
 */

class BasicOrientation {
    private:
        coords local_x; /**< the x, y, z unit vectors for the beam element local x axis. */
        real length; /**< length of the beam-column element. */
        mat T = make_xd_mat(6,12); /**< The 6x12 T transformation matrix. */
        real alpha = 0.0; /**< angle between local and global coordinate systems. */
        real offset = 0.0; /**< offset of centroid along y axis. */
    public:
        /**
         * @brief sets offset and calculates length, local_x axis, alpha, and T matrix
         * 
         * @param nodes a size 2 array of shared ptr to the element nodes
         * @param sec_offset the offset value for the section
         * @param origin_x the global coord system x-axis unit vector
         */
        void evaluate(std::array<std::shared_ptr<Node>, 2> const & nodes, real sec_offset, coords const & origin_x)
        {
            offset = sec_offset;
            calc_length_local_x(nodes);
            calc_alpha(origin_x);
            calc_T();
        }
        /**
         * @brief calculates beam length and local_x components
         * 
         * @param nodes a size 2 array of shared ptr to the element nodes
         */
        void calc_length_local_x(std::array<std::shared_ptr<Node>, 2> const &  nodes) {
            local_x = (nodes[1]->get_coords() - nodes[0]->get_coords());
            length = local_x.norm();
            local_x /= length;
        }
        /**
         * @brief calcualte the angle between the local and global x axes
         * 
         * @details
         * The calculation is based on:
         * \f$ A\dot B = \abs{A}\abs{B}\cos{\theta}\f$
         * and then:
         * \f$ \theta = \arccos{\frac{A\dot B}{\abs{A}\abs{B}}} \f$
         * 
         * @param origin_x the global coord system x-axis unit vector
         */
        void calc_alpha(coords const& origin_x) {
           alpha = std::acos(origin_x.dot(local_x));
        }
        /**
         * @brief 
         * creates the T matrix
         * 
         * @details
         * the T matrix created is a 6x12 matrix with members only in the first
         * 6 rows.
         */
        void calc_T() {
            real c = std::cos(alpha);
            real s = std::sin(alpha);
            T(0,0) = c;
            T(0,2) = s;
            T(0,5) = offset;
            T(1,0) = -s;
            T(1,2) = c;
            T(2,5) = 1;
            T(3,6) = c;
            T(3,8) = s;
            T(3,11) = offset;
            T(4,6) = -s;
            T(4,8) = c;
            T(5,11) = 1;
        }
        mat get_T() {return T;}
};
#endif