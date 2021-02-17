/**
 * @file numpdesetup_test.cc
 * @brief NPDE homework NumPDESetup code
 * @author Oliver Rietmann, Erick Schulz
 * @date 01.01.2020
 * @copyright Developed at ETH Zurich
 */

#include <Eigen/Core>

#include <gtest/gtest.h>

#include "../numpdesetup.h"

namespace NumPDESetup::test {

TEST(NumPDESetup, dummyFunction) {
  double x = 1.0;
  int n = 2;

  Eigen::Vector2d v = NumPDESetup::dummyFunction(x, n);

  Eigen::Vector2d v_ref = {1.0, 1.0};

  double tol = 1.0e-8;
  ASSERT_NEAR(0.0, (v - v_ref).lpNorm<Eigen::Infinity>(), tol);
}

}  // namespace NumPDESetup::test
