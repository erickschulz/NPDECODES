/**
 * @file numpdesetup_test.cc
 * @brief NPDE homework NumPDESetup code
 * @author Ralf Hiptmair, Oliver Rietmann, Erick Schulz
 * @date 17.02.2020
 * @copyright Developed at ETH Zurich
 */

#include "../numpdesetup.h"

#include <gtest/gtest.h>

#include <Eigen/Core>

namespace NumPDESetup::test {

TEST(NumPDESetup, dummyFunction) {
  double x = 1.0;
  int n = 2;

  Eigen::Vector2d v = NumPDESetup::dummyFunction(x, n);
  // The expected solution
  Eigen::Vector2d v_ref = {1.0, 1.0};
  // Check whether function call produces correct result using the Google test
  // framework
  double tol = 1.0e-8;
  ASSERT_NEAR(0.0, (v - v_ref).lpNorm<Eigen::Infinity>(), tol);
}

}  // namespace NumPDESetup::test
