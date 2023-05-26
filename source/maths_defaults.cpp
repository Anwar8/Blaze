#include <cmath>
#include "maths_defaults.hpp"
#ifdef EIGEN
vec make_xd_vec(int rows) {
    return vec::Zero(rows);
}
mat make_xd_mat(int rows, int columns) {
    return mat::Zero(rows, columns);
}
#endif