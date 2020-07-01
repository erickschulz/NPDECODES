/**
 * @file parametricfiniteelements_test.cc
 * @brief NPDE homework ParametricFiniteElements code
 * @author Amélie Loher
 * @date 06.04.2020
 * @copyright Developed at ETH Zurich
 */

#include "../parametricfiniteelements.h"

#include <gtest/gtest.h>

#include <Eigen/Core>
#include <Eigen/Sparse>

#include <vector>

namespace ParametricFiniteElements::test {

TEST(ParametricFiniteElements, assembleGeoTherm) {
  unsigned int n = 2;

  auto alpha = [](Eigen::Vector2d x) -> double { return 3.0 / 2.0; };

  auto Psi = [](double x1) -> double { return x1 * x1 + 1.0; };

  std::vector<Eigen::Triplet<double>> triplets =
      ParametricFiniteElements::assembleGeoTherm(n, alpha, Psi);

  Eigen::SparseMatrix<double> A_sps((n + 1) * (n + 1), (n + 1) * (n + 1));
  A_sps.setFromTriplets(triplets.begin(), triplets.end());

  Eigen::MatrixXd A(A_sps);

  Eigen::MatrixXd A_ref(9, 9);
  A_ref << 1.61719, -0.84375, 0, -0.867188, 0.09375, 0, 0, 0, 0, -0.84375, 3.45,
      -1.21875, -0.09375, -1.575, 0.28125, 0, 0, 0, 0, -1.21875, 1.69922, 0,
      -0.28125, -0.199219, 0, 0, 0, -0.867188, -0.09375, 0, 3.32812, -1.6875, 0,
      -0.960938, 0.28125, 0, 0.09375, -1.575, -0.28125, -1.6875, 7.65, -2.4375,
      -0.28125, -2.325, 0.84375, 0, 0.28125, -0.199219, 0, -2.4375, 3.82031, 0,
      -0.84375, -0.621094, 0, 0, 0, -0.960938, -0.28125, 0, 2.36545, -0.84375,
      0, 0, 0, 0, 0.28125, -2.325, -0.84375, -0.84375, 5.6802, -1.21875, 0, 0,
      0, 0, 0.84375, -0.621094, 0, -1.21875, 1.44679;

  double tol = 1.0e-4;

  ASSERT_NEAR(0.0, (A - A_ref).lpNorm<Eigen::Infinity>(), tol);
}

TEST(ParametricFiniteElements, geoThermSolve) {
  unsigned int n = 3;

  auto alpha = [](Eigen::Vector2d x) -> double { return 3.0 / 2.0; };

  auto Psi = [](double x1) -> double { return x1 * x1 + 1.0; };

  Eigen::VectorXd mu = ParametricFiniteElements::geoThermSolve(n, alpha, Psi);

  Eigen::VectorXd mu_ref((n + 1) * (n + 1));
  mu_ref << 1, 1, 1, 1, 0.860946, 0.845351, 0.79932, 0.724646, 0.722159,
      0.690899, 0.601355, 0.460807, 0.584267, 0.537488, 0.411906, 0.220445;

  double tol = 1.0e-6;

  ASSERT_NEAR(0.0, (mu - mu_ref).lpNorm<Eigen::Infinity>(), tol);
}

TEST(ParametricFiniteElements, geoThermSurfInt) {
  unsigned int n = 3;

  auto Psi = [](double x1) -> double { return x1 * x1 + 1.0; };

  Eigen::VectorXd mu_ref((n + 1) * (n + 1));
  mu_ref << 1, 1, 1, 1, 0.860946, 0.845351, 0.79932, 0.724646, 0.722159,
      0.690899, 0.601355, 0.460807, 0.584267, 0.537488, 0.411906, 0.220445;

  double val = geoThermSurfInt(n, Psi, mu_ref);

  double val_ref = 0.625691;

  double tol = 1.0e-6;
  ASSERT_NEAR(0.0, val - val_ref, tol);
}

}  // namespace ParametricFiniteElements::test
