/**
 * @file crossprod_test.cc
 * @brief NPDE homework CrossProd code
 * @author Oliver Rietmann
 * @date 31.03.2021
 * @copyright Developed at ETH Zurich
 */

#include "../crossprod.h"

#include <gtest/gtest.h>

#include <Eigen/Core>
#include <iostream>

namespace CrossProd::test {

Eigen::Vector3d f(Eigen::Vector3d y) {
  return Eigen::Vector3d(y(0) * y(1), y(1) * y(2), y(2) - y(0));
};

Eigen::Matrix3d Jf(Eigen::Vector3d y) {
  Eigen::Matrix3d J;
  J << y(1), y(0), 0.0, 0.0, y(2), y(1), -1.0, 0.0, 1.0;
  return J;
};

TEST(CrossProd, solve_imp_mid) {
  double T = 1.0;
  int N = 1;
  Eigen::Vector3d y0(0.1, 0.2, 0.4);
  Eigen::Vector3d yT = solve_imp_mid(f, Jf, T, y0, N).back();

  Eigen::Vector3d yT_reference(0.135781669194583, 0.407030550935246,
                               0.964218330194744);

  double tol = 1.0e-5;
  ASSERT_NEAR(0.0, (yT - yT_reference).lpNorm<Eigen::Infinity>(), tol);
}

TEST(CrossProd, solve_lin_mid) {
  double T = 1.0;
  int N = 1;
  Eigen::Vector3d y0(0.1, 0.2, 0.4);
  Eigen::Vector3d yT = solve_lin_mid(f, Jf, T, y0, N).back();

  Eigen::Vector3d yT_reference(0.131724137931034, 0.371034482758621,
                               0.968275862068966);

  double tol = 1.0e-8;
  ASSERT_NEAR(0.0, (yT - yT_reference).lpNorm<Eigen::Infinity>(), tol);
}

}  // namespace CrossProd::test
