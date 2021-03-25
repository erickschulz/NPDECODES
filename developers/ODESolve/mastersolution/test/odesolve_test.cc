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

double Psi(double h, double y0) { return y0 * (1 + h); }

TEST(ODESolve, PsiTilde) {
  double y0 = 1.0;
  double result = PsiTilde(Psi, 1, 0.1, y0);
  double reference = 1.105;
  double tol = 1.0e-8;
  ASSERT_NEAR(reference, result, tol);
}

TEST(ODESolve, OdeIntEqui) {
  double y0 = 1.0;
  double T = 1.0;
  int M = 8;

  Eigen::VectorXd result(M + 1);
  std::vector<double> vector = OdeIntEqui(Psi, T, y0, M);
  for (int m = 0; m < M + 1; ++m) {
    result(m) = vector[m];
  }

  Eigen::VectorXd reference(M + 1);
  reference << 1.0, 1.125, 1.265625, 1.423828125, 1.601806640625,
      1.80203247070312, 2.02728652954102, 2.28069734573364, 2.56578451395035;

  double tol = 1.0e-8;
  Eigen::VectorXd error = reference - result;
  ASSERT_NEAR(0.0, error.lpNorm<Eigen::Infinity>(), tol);
}

TEST(ODESolve, OdeIntSsCtrl) {
  double y0 = 1.0;
  double T = 1.0;
  auto [dummy, vector] = OdeIntSsCtrl(Psi, 1, y0, 1, 0.01, 10e-5, 10e-5, 10e-5);

  Eigen::Vector2d result(vector[0], vector[1]);
  Eigen::Vector2d reference(1.0, 1.01005);

  double tol = 1.0e-8;
  Eigen::Vector2d error = reference - result;
  ASSERT_NEAR(0.0, error.lpNorm<Eigen::Infinity>(), tol);
}

}  // namespace ODESolve::test
