#ifndef MATHS_DEFAULTS_HPP
#define MATHS_DEFAULTS_HPP

#include <iostream>

using real = double;
using realx2 = long double;

#define EIGEN

#ifdef EIGEN
#include <Eigen/Dense>
using coords = Eigen::Vector3d;
using vec = Eigen::VectorXd;
using mat = Eigen::MatrixXd;
#endif

vec make_xd_vec(int rows);
mat make_xd_mat(int rows, int columns);

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
