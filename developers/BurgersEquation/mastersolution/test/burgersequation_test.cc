/**
 * @file burgersequation_test_mastersolution.cc
 * @brief NPDE homework BurgersEquation code
 * @author Oliver Rietmann
 * @date 15.04.2019
 * @copyright Developed at ETH Zurich
 */

#include "../burgersequation.h"

#include <gtest/gtest.h>

#include <Eigen/Core>

namespace BurgersEquation::test {

TEST(BurgersEquation, solveBurgersGodunov) {
  double T = 2.0;
  unsigned int N = 10;

  Eigen::VectorXd mu = BurgersEquation::solveBurgersGodunov(T, N);

  Eigen::VectorXd mu_ref(N + 1);
  mu_ref << 0, 0, 0, 0.10703626063314, 0.228642274632239, 0.345805344695093,
      0.281425975977912, 0.0370901440616156, 0, 0, 0;

  double tol = 1.0e-8;
  ASSERT_NEAR(0.0, (mu - mu_ref).lpNorm<Eigen::Infinity>(), tol);
}

TEST(BurgersEquation, numexpBurgersGodunov) {
  Eigen::Matrix<double, 3, 4> result = BurgersEquation::numexpBurgersGodunov();

  Eigen::Matrix<double, 3, 4> result_ref;
  result_ref << 0.1, 0.05, 0.025, 0.0125, 0.0530488551171146,
      0.0263066209538425, 0.0120042822755487, 0.00584742030018541,
      0.0644629151855385, 0.0401571668250919, 0.0183165195121137,
      0.00779107247088863;

  double tol = 1.0e-8;
  ASSERT_NEAR(0.0, (result - result_ref).lpNorm<Eigen::Infinity>(), tol);
}

}  // namespace BurgersEquation::test
