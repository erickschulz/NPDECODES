/**
 * @file gradientflow_test.cc
 * @brief NPDE homework GradientFlow code
 * @author Oliver Rietmann
 * @date 04.04.2021
 * @copyright Developed at ETH Zurich
 */

#include "../gradientflow.h"

#include <gtest/gtest.h>

#include <Eigen/Core>
#include <array>
#include <vector>

namespace GradientFlow::test {

Eigen::Vector2d f(Eigen::Vector2d x) {
  return Eigen::Vector2d(x(0) * x(0) + x(1) * x(1), x(0) * x(1));
}

Eigen::Matrix2d df(Eigen::Vector2d x) {
  Eigen::Matrix2d J;
  J << 2 * x(0), 2 * x(1), x(1), x(0);
  return J;
}

TEST(GradientFlow, computeStages) {
  double h = 0.5;
  Eigen::Vector2d y0(1.0, 0.5);

  std::array<Eigen::VectorXd, 5> stages = ComputeStages(f, df, y0, h);
  Eigen::MatrixXd G(5, 2);
  for (int i = 0; i < stages.size(); ++i) {
    G.row(i) = stages[i];
  }

  Eigen::MatrixXd G_reference(5, 2);
  G_reference << 0.234048464188, 0.0911954850442, 1.16124077723, 0.43498079856,
      0.589539233307, 0.228323917228, 0.499164120689, 0.194780654616,
      1.04447583358, 0.415950635186;

  double tol = 1.0e-6;
  ASSERT_NEAR(0.0, (G - G_reference).lpNorm<Eigen::Infinity>(), tol);
}

TEST(GradientFlow, DiscEvolSDIRK) {
  double h = 0.5;
  Eigen::Vector2d y0(1.0, 0.5);

  Eigen::VectorXd yh = DiscEvolSDIRK(f, df, y0, h);
  ASSERT_EQ(yh.size(), 2);

  Eigen::Vector2d yh_reference(2.04447583358, 0.915950635186);

  double tol = 1.0e-6;
  ASSERT_NEAR(0.0, (yh - yh_reference).lpNorm<Eigen::Infinity>(), tol);
}

constexpr double SQRT2 = 1.41421356237309504880;

TEST(GradientFlow, SolveGradientFlow) {
  Eigen::Vector2d d(0.5 * SQRT2, 0.5 * SQRT2);
  double lambda = 8.5;
  Eigen::Vector2d y0(1.0, 0.0);
  double T = 0.5;
  int M = 10;

  std::vector<Eigen::VectorXd> yvector = SolveGradientFlow(d, lambda, y0, T, M);
  Eigen::MatrixXd Y(M + 1, 2);
  for (int i = 0; i < yvector.size(); ++i) {
    Y.row(i) = yvector[i];
  }

  Eigen::MatrixXd Y_reference(M + 1, 2);
  Y_reference << 1.0, 0.0, 0.661574350863, -0.265183334169, 0.500522309998,
      -0.345715652452, 0.415074659032, -0.354835307771, 0.36129317252,
      -0.337895066294, 0.321645925272, -0.312568326946, 0.289168885473,
      -0.28564991388, 0.26103327209, -0.259669864129, 0.236018396524,
      -0.235490344277, 0.213530281589, -0.21332581582, 0.193223336322,
      -0.193144178833;

  double tol = 1.0e-6;
  ASSERT_NEAR(0.0, (Y - Y_reference).lpNorm<Eigen::Infinity>(), tol);
}

}  // namespace GradientFlow::test
