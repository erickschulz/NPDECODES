/**
 * @file exponentialintegrator.cc
 * @brief NPDE homework ExponentialIntegrator code
 * @author Unknown, Oliver Rietmann
 * @date 04.04.2021
 * @copyright Developed at ETH Zurich
 */

#include "exponentialintegrator.h"

#include <Eigen/Core>
#include <cassert>
#include <unsupported/Eigen/MatrixFunctions>

namespace ExponentialIntegrator {

//! \brief Function $\phi$ used in the Exponential Euler
//! single step method for an autonomous ODE.
Eigen::MatrixXd phim(const Eigen::MatrixXd &Z) {
  int n = Z.cols();
  assert(n == Z.rows() && "Matrix must be square.");
  Eigen::MatrixXd C(2 * n, 2 * n);
  C << Z, Eigen::MatrixXd::Identity(n, n), Eigen::MatrixXd::Zero(n, 2 * n);
  return C.exp().block(0, n, n, n);
}

}  // namespace ExponentialIntegrator
