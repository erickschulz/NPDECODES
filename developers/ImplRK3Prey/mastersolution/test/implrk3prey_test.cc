/**
 * @file implrk3prey_test.cc
 * @brief NPDE homework ImplRK3Prey code
 * @author Oliver Rietmann
 * @date 29.03.2021
 * @copyright Developed at ETH Zurich
 */

#include "../implrk3prey.h"

#include <gtest/gtest.h>

#include <Eigen/Core>

namespace ImplRK3Prey::test {

TEST(ImplRK3Prey, solve) {
  // Set up ODE for predator/prey model
  Eigen::Vector2d alpha(3.0, 2.0);
  Eigen::Vector2d beta(0.1, 0.1);
  auto f = [alpha, beta] (Eigen::Vector2d y) -> Eigen::Vector2d {
    return Eigen::Vector2d(y(0) * (alpha(0) - beta(0) * y(1)), y(1) * (-alpha(1) + beta(1) * y(0)));
  };
  auto Jf = [alpha, beta] (Eigen::Vector2d y) -> Eigen::Matrix2d {
	  Eigen::Matrix2d J;
    J << alpha(0) - beta(0) * y(1), -beta(0) * y(0), beta(1) * y(1), -alpha(1) + beta(1) * y(0);
    return J;
  };

  // Set up the integrator
  Eigen::Matrix2d A;
  A << 5.0 / 12.0, -1.0 / 12.0, 3.0 / 4.0, 1.0 / 4.0;
  Eigen::Vector2d b(0.75, 0.25);
  implicitRKIntegrator integrator(A, b);

  // Perform the integration, i.e. compute y(T)
  double T = 2.0;
  int N = 128;
  Eigen::Vector2d y0(100.0, 5.0);
  Eigen::Vector2d yT = integrator.solve(f, Jf, T, y0, N).back();

  // Compare with reference solution
  Eigen::Vector2d yT_reference(0.372657928550525, 8.41720603055351);
  double tol = 1.0e-8;
  ASSERT_NEAR(0.0, (yT_reference - yT).lpNorm<Eigen::Infinity>(), tol);
}

}  // namespace ImplRK3Prey::test
