/**
 * @file sdirk_test.cc
 * @brief NPDE homework SDIRK code
 * @author Oliver Rietmann
 * @date 31.03.2021
 * @copyright Developed at ETH Zurich
 */

#include "../sdirk.h"

#include <gtest/gtest.h>

#include <Eigen/Core>
#include <cmath>
#include <iostream>

namespace SDIRK::test {

TEST(SDIRK, sdirkStep) {
  double T = 10.0;
  int N = 320;
  double gamma = (3.0 + std::sqrt(3.0)) / 6.0;
  Eigen::Vector2d z0(1.0, 0.0);

  double h = T / N;
  Eigen::Vector2d zT = SdirkStep(z0, h, gamma);

  Eigen::Vector2d zT_reference(0.999516810128635, -0.0307616782094714);

  double tol = 1.0e-8;
  ASSERT_NEAR(0.0, (zT - zT_reference).lpNorm<Eigen::Infinity>(), tol);
}

TEST(SDIRK, sdirkSolve) {
  double T = 10.0;
  int N = 320;
  double gamma = (3.0 + std::sqrt(3.0)) / 6.0;
  Eigen::Vector2d z0(1.0, 0.0);

  Eigen::Vector2d zT = SdirkSolve(z0, N, T, gamma).back();

  Eigen::Vector2d zT_reference(-0.00216997617476873, -0.0053856854265896);

  double tol = 1.0e-8;
  ASSERT_NEAR(0.0, (zT - zT_reference).lpNorm<Eigen::Infinity>(), tol);
}

}  // namespace SDIRK::test
