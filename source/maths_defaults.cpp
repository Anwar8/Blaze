#include <cmath>
#include "maths_defaults.hpp"
#ifdef EIGEN
vec make_xd_vec(int rows) {
    return vec::Zero(rows);
}
mat make_xd_mat(int rows, int columns) {
    return mat::Zero(rows, columns);
}

spvec make_spd_vec(int rows) {
    return spvec(rows);
}
spmat make_spd_mat(int rows, int columns) {
    return spmat(rows, columns);
}
#endif
