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

TEST(NonLinSchroedingerEquation, MassElementMatrixProvider) {
  Eigen::Matrix3d D = Eigen::Matrix3d::Identity();
  Eigen::Matrix3d D_ref = Eigen::Matrix3d::Identity();
  double tol = 1.0e-8;
  ASSERT_NEAR(0.0, (D - D_ref).lpNorm<Eigen::Infinity>(), tol);
}

}  // namespace NonLinSchroedingerEquation::test
