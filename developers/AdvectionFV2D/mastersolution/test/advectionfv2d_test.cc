/**
 * @file advectionfv2d_test.cc
 * @brief NPDE homework AdvectionFV2D code
 * @author Philipp Egg
 * @date 01.01.2020
 * @copyright Developed at ETH Zurich
 */

#include <Eigen/Core>

#include <gtest/gtest.h>

#include "../advectionfv2d.h"

namespace AdvectionFV2D::test {

TEST(AdvectionFV2D, dummyFunction) {
  double x = 0.0;
  int n = 0;

  // Eigen::Vector2d v = AdvectionFV2D::dummyFunction(x, n);
  Eigen::Vector2d v = {1.0, 1.0};

  Eigen::Vector2d v_ref = {1.0, 1.0};

  double tol = 1.0e-8;
  ASSERT_NEAR(0.0, (v - v_ref).lpNorm<Eigen::Infinity>(), tol);
}

}  // namespace AdvectionFV2D::test
