/**
 * @file  maximumprinciple_test.cc
 * @brief NPDE homework "MaximumPrinciple" code
 * @author Oliver Rietmann
 * @date 25.03.2019
 * @copyright Developed at ETH Zurich
 */

#include "../maximumprinciple.h"

#include <gtest/gtest.h>

#include <Eigen/Core>
#include <functional>

namespace MaximumPrinciple::test {

TEST(MaximumPrinciple, computeGalerkinMatrix) {
  unsigned int M = 3;
  int M2 = M * M;
  double c;
  Eigen::MatrixXd A(M2, M2);
  Eigen::MatrixXd A_ref(M2, M2);

  double tol = 1e-8;

  c = 0.0;
  A = Eigen::MatrixXd(computeGalerkinMatrix(M, c));
  A_ref << 4, -1, 0, -1, 0, 0, 0, 0, 0, -1, 4, -1, 0, -1, 0, 0, 0, 0, 0, -1, 4,
      0, 0, -1, 0, 0, 0, -1, 0, 0, 4, -1, 0, -1, 0, 0, 0, -1, 0, -1, 4, -1, 0,
      -1, 0, 0, 0, -1, 0, -1, 4, 0, 0, -1, 0, 0, 0, -1, 0, 0, 4, -1, 0, 0, 0, 0,
      0, -1, 0, -1, 4, -1, 0, 0, 0, 0, 0, -1, 0, -1, 4;
  ASSERT_NEAR(0.0, (A - A_ref).lpNorm<Eigen::Infinity>(), tol);

  c = 1.0;
  A = Eigen::MatrixXd(computeGalerkinMatrix(M, c));
  A_ref << 0.03125, 0.00520833, 0, 0.00520833, 0.00520833, 0, 0, 0, 0,
      0.00520833, 0.03125, 0.00520833, 0, 0.00520833, 0.00520833, 0, 0, 0, 0,
      0.00520833, 0.03125, 0, 0, 0.00520833, 0, 0, 0, 0.00520833, 0, 0, 0.03125,
      0.00520833, 0, 0.00520833, 0.00520833, 0, 0.00520833, 0.00520833, 0,
      0.00520833, 0.03125, 0.00520833, 0, 0.00520833, 0.00520833, 0, 0.00520833,
      0.00520833, 0, 0.00520833, 0.03125, 0, 0, 0.00520833, 0, 0, 0, 0.00520833,
      0, 0, 0.03125, 0.00520833, 0, 0, 0, 0, 0.00520833, 0.00520833, 0,
      0.00520833, 0.03125, 0.00520833, 0, 0, 0, 0, 0.00520833, 0.00520833, 0,
      0.00520833, 0.03125;
  ASSERT_NEAR(0.0, (A - A_ref).lpNorm<Eigen::Infinity>(), tol);
}

TEST(MaximumPrinciple, computeLoadVector) {
  unsigned int M = 3;
  int M2 = M * M;

  double tol = 1e-6;

  auto f = [](double x, double y) { return -(x * x + y * y); };

  Eigen::VectorXd phi = computeLoadVector(M, f);
  Eigen::VectorXd phi_ref(M2);
  phi_ref << -0.0078125, -0.01953125, -0.0390625, -0.01953125, -0.03125,
      -0.05078125, -0.0390625, -0.05078125, -0.0703125;

  for (int i = 0; i < M2; ++i) {
    EXPECT_NEAR(phi(i), phi_ref(i), tol);
  }
}

TEST(MaximumPrinciple, computeGalerkinMatrixTR) {
  unsigned int M = 3;
  int M2 = M * M;
  double c;
  Eigen::MatrixXd A(M2, M2);
  Eigen::MatrixXd A_ref(M2, M2);

  double tol = 1e-8;

  c = 0.0;
  A = Eigen::MatrixXd(computeGalerkinMatrixTR(M, c));
  A_ref << 4, -1, 0, -1, 0, 0, 0, 0, 0, -1, 4, -1, 0, -1, 0, 0, 0, 0, 0, -1, 4,
      0, 0, -1, 0, 0, 0, -1, 0, 0, 4, -1, 0, -1, 0, 0, 0, -1, 0, -1, 4, -1, 0,
      -1, 0, 0, 0, -1, 0, -1, 4, 0, 0, -1, 0, 0, 0, -1, 0, 0, 4, -1, 0, 0, 0, 0,
      0, -1, 0, -1, 4, -1, 0, 0, 0, 0, 0, -1, 0, -1, 4;
  ASSERT_NEAR(0.0, (A - A_ref).lpNorm<Eigen::Infinity>(), tol);

  c = 1.0;
  A = Eigen::MatrixXd(computeGalerkinMatrixTR(M, c));
  A_ref = 0.0625 * Eigen::MatrixXd::Identity(M2, M2);
  ASSERT_NEAR(0.0, (A - A_ref).lpNorm<Eigen::Infinity>(), tol);
}

}  // namespace MaximumPrinciple::test
