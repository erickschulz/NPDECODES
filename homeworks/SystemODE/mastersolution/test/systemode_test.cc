/**
 * @file systemode_test.cc
 * @brief NPDE homework SystemODE
 * @copyright Developed at ETH Zurich
 */

#include "../systemode.h"

#include <gtest/gtest.h>

#include <Eigen/Core>

namespace SystemODE::test {

TEST(SystemODE, rk4step) {
  double h = 0.5;  // stepsize

  // Right-hand side of ode
  auto f = [](Eigen::VectorXd y) {
    Eigen::VectorXd eval(3);
    eval << y(1) * y(2), y(0) * y(1), 3.0 * y(2);
    return eval;
  };

  // Initial value
  Eigen::VectorXd y0(3);
  y0 << -1.0, 1.0, 2.0;

  // RK4 step
  Eigen::VectorXd result = SystemODE::rk4step(f, h, y0);

  Eigen::Vector3d reference(0.931516011555989, 0.879332304000854, 8.796875);

  double tol = 1.0e-8;
  Eigen::Vector3d error = reference - result;
  ASSERT_NEAR(0.0, error.lpNorm<Eigen::Infinity>(), tol);
}

}  // namespace SystemODE::test
