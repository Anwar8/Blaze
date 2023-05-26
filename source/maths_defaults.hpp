#ifndef MATHS_DEFAULTS_HPP
#define MATHS_DEFAULTS_HPP



using real = double;
using realx2 = long double;

#define EIGEN

#ifdef EIGEN
#include <Eigen/Dense>
using coords = Eigen::Vector3d;
using vec = Eigen::VectorXd;
using mat = Eigen::MatrixXd;
#endif


#endif
