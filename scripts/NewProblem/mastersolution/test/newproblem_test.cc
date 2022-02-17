/**
 * @file newproblem_test.cc
 * @brief NPDE homework NewProblem code
 * @author Oliver Rietmann, Erick Schulz
 * @date 01.01.2020
 * @copyright Developed at ETH Zurich
 */

#include "../newproblem.h"

#include <gtest/gtest.h>

#include <Eigen/Core>

namespace NewProblem::test {

TEST(NewProblem, dummyFunction) {
  double x = 0.0;
  int n = 0;

  Eigen::Vector2d v = NewProblem::dummyFunction(x, n);

  Eigen::Vector2d v_ref = {1.0, 1.0};

  double tol = 1.0e-8;
  ASSERT_NEAR(0.0, (v - v_ref).lpNorm<Eigen::Infinity>(), tol);
}

}  // namespace NewProblem::test
