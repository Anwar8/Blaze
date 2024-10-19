/**
 * @file maths_defaults.hpp
 * @brief Contains wrappers for numerical types, containers, and operations. 
 * 
 */
#ifndef MATHS_DEFAULTS_HPP
#define MATHS_DEFAULTS_HPP

#include <iostream>
/**
 * @defgroup MathsTypesFunctions
 * @brief The types and functions used for representing mathematical types and operations.
 * @{
 */
using real = double; /**< Numerical data type to represent double precision real numbers.*/
using realx2 = long double; /**< Numerical data type to represent quadruple precision real numbers.*/
/**
 * @brief Numerical engine being used.
 * 
 * @details Allowed options are: EIGEN.
 * 
 */
#define EIGEN /**< Numerical engine being used.*/

#ifdef EIGEN
#include <Eigen/Dense>
#include <Eigen/Sparse>

/**
 * @brief Alias for a \ref real vector with three rows. 
 */
using coords = Eigen::Matrix<real, 3, 1>;
/**
 * @brief Alias for a dense \ref real vector with dynamic rows. 
 */
using vec = Eigen::Matrix<real, Eigen::Dynamic, 1>;
/**
 * @brief Alias for a dense \ref real matrix with dynamic rows and columns. 
 */
using mat = Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic>;

// sparse non-zero element: val, row, col
/**
 * @brief Alias for triplet-indexed \ref real numbers to be used for sparse matrix creation.
 */
using spnz = Eigen::Triplet<real>;
/**
 * @brief Alias for a sparse \ref real vector.
 */
using spvec = Eigen::SparseVector<real>;
/**
 * @brief Alias for a sparse \ref real matrix.
 */
using spmat = Eigen::SparseMatrix<real>;
#endif
/**
 * @brief Allocates a dense vector with dynamically-allocated rows.
 * 
 * @param rows Number of rows in vector.
 * @return vec Dense vector allocated.
 */
inline vec make_xd_vec(int rows)
{
    return vec::Zero(rows);
}
/**
 * @brief Allocates a dense matrix with dynamically-allocated rows and columns.
 * 
 * @param rows Number of rows in matrix.
 * @param columns Number of columns in matrix.
 * @return mat Dense matrix allocated.
 */
inline mat make_xd_mat(int rows, int columns)
{
    return mat::Zero(rows, columns);
}
/**
 * @brief Allocates a sparse vector with given number of rows.
 * 
 * @param rows Number of rows to allocate to sparse vector.
 * @return spvec Sparse vector allocated.
 */
inline spvec make_spd_vec(int rows)
{
    return spvec(rows);
}
/**
 * @brief Allocates a sparse matrix with given number of rows and columns.
 * 
 * @param rows Number of rows to allocate to sparse matrix.
 * @param columns Number of columns to allocate to sparse matrix.
 * @return spmat Sparse matrix allocated.
 */
inline spmat make_spd_mat(int rows, int columns)
{
    return spmat(rows, columns);
}
/**
 * @brief return the l2 norm of a vector.
 * 
 * @tparam eigen_container any container form Eigen compatible with `squaredNorm()`.
 * @param V the vector/container to be normalised.
 * @return realx2 l2_norm twice the normal precision norm value.
 */
template <typename eigen_container>
inline realx2 calc_l2_norm(eigen_container& V)
{
    realx2 l2_norm = V.squaredNorm();
    return l2_norm;
}
/** @} */

#endif
