/**
 * @file nonlinschroedingerequation_test.cc
 * @brief NPDE homework NonLinSchroedingerEquation code
 * @author Oliver Rietmann
 * @date 22.04.2020
 * @copyright Developed at ETH Zurich
 */

#include <Eigen/Core>

#include <gtest/gtest.h>

#include "../nonlinschroedingerequation.h"

namespace NonLinSchroedingerEquation::test {

TEST(NonLinSchroedingerEquation, dummyFunction) {
  double x = 0.0;
  int n = 0;

  Eigen::Vector2d v = NonLinSchroedingerEquation::dummyFunction(x, n);

  Eigen::Vector2d v_ref = {1.0, 1.0};

  double tol = 1.0e-8;
  ASSERT_NEAR(0.0, (v - v_ref).lpNorm<Eigen::Infinity>(), tol);
}

}  // namespace NonLinSchroedingerEquation::test
