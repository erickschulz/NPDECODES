/**
 * @file rk3prey_test.cc
 * @brief NPDE homework RK3Prey code
 * @author Oliver Rietmann
 * @date 14.03.2021
 * @copyright Developed at ETH Zurich
 */

#include "../rkintegrator.h"

#include <vector>

#include <gtest/gtest.h>

#include <Eigen/Core>

namespace RK3Prey::test {

TEST(RK3Prey, RKIntegrator) {
  Eigen::Matrix2d A;
  A << 0, 0, 1, 0;
  Eigen::Vector2d b(0.5, 0.5);
  RKIntegrator<Eigen::VectorXd> rkintegrator(A, b);

	double T = 2.0;
  int N = 10;
  Eigen::Vector2d y0(-1.0, 1.0);
  auto f = [] (Eigen::VectorXd y) { return Eigen::Vector2d(-0.5 * y(0), y(0) * y(1)); };
  std::vector<Eigen::VectorXd> result = rkintegrator.solve(f, T, y0, N);

  std::vector<Eigen::VectorXd> reference(N + 1);
  reference[0] = Eigen::Vector2d(-1, 1);
  reference[1] = Eigen::Vector2d(-0.905, 0.828);
  reference[2] = Eigen::Vector2d(-0.819025, 0.6978321486);
  reference[3] = Eigen::Vector2d(-0.741217625, 0.597665102250463);
  reference[4] = Eigen::Vector2d(-0.670801950625, 0.519405587909154);
  reference[5] = Eigen::Vector2d(-0.607075765315625, 0.457413068349474);
  reference[6] = Eigen::Vector2d(-0.54940356761064, 0.40768739285642);
  reference[7] = Eigen::Vector2d(-0.49721022868763, 0.367345306253303);
  reference[8] = Eigen::Vector2d(-0.449975256962305, 0.334276874715393);
  reference[9] = Eigen::Vector2d(-0.407227607550886, 0.306916078643447);
  reference[10] = Eigen::Vector2d(-0.368540984833552, 0.284085135532343);

  Eigen::VectorXd error(N + 1);
  for (int n = 0; n < N + 1; ++n) {
    error(n) = (reference[n] - result[n]).lpNorm<Eigen::Infinity>();
  }
  double tol = 1.0e-8;
  ASSERT_NEAR(0.0, error.lpNorm<Eigen::Infinity>(), tol);
}

}  // namespace RK3Prey::test
