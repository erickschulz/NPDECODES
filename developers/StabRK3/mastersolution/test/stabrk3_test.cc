/**
 * @file stabrk3_test.cc
 * @brief NPDE homework StabRK3 code
 * @author Philippe Peter
 * @date 13.04.2021
 * @copyright Developed at ETH Zurich
 */

#include "../stabrk3.h"

#include <gtest/gtest.h>

#include <Eigen/Core>

namespace StabRK3::test {

TEST(StabRK3, PredPrey) {
  double T = 1.0;
  Eigen::Vector2d y0(100.0, 1.0);
  double N = 128;
  Eigen::Vector2d yN = PredPrey(y0, T, N);

  // reference solution:
  Eigen::Vector2d y_ref(6.47562991232002e-27, 40.7385871913941);
  double error = (yN - y_ref).norm();

  EXPECT_NEAR(error, 0.0, 1.0E-5);
}

}  // namespace StabRK3::test