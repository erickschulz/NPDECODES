/**
 * @file systemode_test.cc
 * @brief NPDE homework SystemODE code
 * @author Oliver Rietmann
 * @date 14.03.2021
 * @copyright Developed at ETH Zurich
 */

#include "../systemode.h"

#include <gtest/gtest.h>

#include <Eigen/Core>

namespace SystemODE::test {

Eigen::Vector3d f(Eigen::Vector3d y) {
  return Eigen::Vector3d(y(1) * y(2), y(0) * y(1), 3.0 * y(2));
}

TEST(SystemODE, rk4step) {
  double h = 0.5;
  Eigen::Vector3d y0(-1.0, 1.0, 2.0);
  Eigen::Vector3d result;
  rk4step(f, h, y0, result);

  Eigen::Vector3d reference(0.931516011555989, 0.879332304000854, 8.796875);

  double tol = 1.0e-8;
  Eigen::Vector3d error = reference - result;
  ASSERT_NEAR(0.0, error.lpNorm<Eigen::Infinity>(), tol);
}

}  // namespace SystemODE::test
