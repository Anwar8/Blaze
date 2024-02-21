/**
 * @file numerical_debuggers.hpp
 * @brief debuggers to find numerical issues that prevent convergence
 * 
 */
#ifndef NUMERICAL_DEBUGGERS
#define NUMERICAL_DEBUGGERS

#include "maths_defaults.hpp"

/**
 * @defgroup NumericalDebuggers
 * @brief checks if matrices are ill-conditioned or have other issues that prevent convergence.
 * @{
 */

/**
 * @brief Supposed to check if a matrix has numerical issues such as a zero row or a zero diagonal.
 * 
 * @warning does NOT work.
 * 
 * @todo Fix.
 * 
 * @param A sparse matrix to check.
 * @return true if matrix has a numerical error.
 * @return false if matrix does not have a numerical error.
 */
bool check_matrix(spmat A);

/**
 * @brief check if a sparse matrix has a zero row.
 * 
 * @warning does NOT work.
 * 
 * @todo Fix.
 * 
 * @param A sparse matrix to check.
 * @return true matrix has zero row.
 * @return false matrix does not have a zero row.
 */
bool has_zero_row(spmat A);
/** @} */ // end of NumericalDebuggers group
#endif