#ifndef MATHS_DEFAULTS_HPP
#define MATHS_DEFAULTS_HPP

#include <iostream>

using real = double;
using realx2 = long double;

#define EIGEN

#ifdef EIGEN
#include <Eigen/Dense>
#include <Eigen/Sparse>
using coords = Eigen::Vector3d;
using vec = Eigen::VectorXd;
using mat = Eigen::MatrixXd;

// sparse non-zero element: val, row, col
using spnz = Eigen::Triplet<real>;
using spvec = Eigen::SparseVector<real>;
using spmat = Eigen::SparseMatrix<real>;
#endif

vec make_xd_vec(int rows);
mat make_xd_mat(int rows, int columns);
spvec make_spd_vec(int rows);
spmat make_spd_mat(int rows, int columns);

template <typename T>
void print_container(T V)
{
  for (auto v: V)
  {
    std::cout << v << " ";
  }
  std::cout << std::endl;
}


#endif
