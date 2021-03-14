/**
 * @file odesolve_test.cc
 * @brief NPDE homework ODESolve code
 * @author Oliver Rietmann
 * @date 14.03.2021
 * @copyright Developed at ETH Zurich
 */

#include "../odesolve.h"

#include <gtest/gtest.h>

#include <Eigen/Core>
#include <vector>

namespace ODESolve::test {

Eigen::VectorXd Psi(double h, const Eigen::VectorXd& y0) {
  return y0 * (1 + h);
}

TEST(ODESolve, psitilde) {
  Eigen::VectorXd y0(1);
  y0 << 1.0;

  Eigen::VectorXd result = psitilde(Psi, 1, 0.1, y0);

  double reference = 1.105;

  double tol = 1.0e-8;
  ASSERT_NEAR(reference, result[0], tol);
}

TEST(ODESolve, odeintequi) {
  Eigen::VectorXd y0(1);
  y0 << 1.0;
  double T = 1.0;
  int N = 8;

  Eigen::VectorXd result(N + 1);
  std::vector<Eigen::VectorXd> vector = odeintequi(Psi, T, y0, N);
  for (int n = 0; n < N + 1; ++n) {
    result(n) = vector[n](0);
  }

  Eigen::VectorXd reference(N + 1);
  reference << 1.0, 1.125, 1.265625, 1.423828125, 1.601806640625,
      1.80203247070312, 2.02728652954102, 2.28069734573364, 2.56578451395035;

  double tol = 1.0e-8;
  Eigen::VectorXd error = reference - result;
  ASSERT_NEAR(0.0, error.lpNorm<Eigen::Infinity>(), tol);
}

TEST(ODESolve, odeintssctrl) {
  Eigen::VectorXd y0(1);
  y0 << 1.0;
  double T = 1.0;
  std::vector<Eigen::VectorXd> vector =
      odeintssctrl(Psi, T, y0, 0.01, 1, 10e-5, 10e-5, 10e-5);
  Eigen::Vector2d result(vector[0](0), vector[1](0));

  Eigen::Vector2d reference(1.0, 1.01005);

  double tol = 1.0e-8;
  Eigen::Vector2d error = reference - result;
  ASSERT_NEAR(0.0, error.lpNorm<Eigen::Infinity>(), tol);
}

}  // namespace ODESolve::test
