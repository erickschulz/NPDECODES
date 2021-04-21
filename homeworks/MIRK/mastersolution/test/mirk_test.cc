/**
 * @file mirk_test.cc
 * @brief NPDE homework MIRK code
 * @author Philippe Peter
 * @date 14.04.2021
 * @copyright Developed at ETH Zurich
 */

#include "../mirk.h"

#include <gtest/gtest.h>

#include <Eigen/Core>

namespace MIRK::test {

TEST(MIRK, Newton2Steps) {
  Eigen::Vector2d x(-1.0, 1.0);

  // f(x,y) = [e^x,e^y]
  auto f = [](Eigen::Vector2d x) {
    return Eigen::Vector2d(std::exp(x(0)), std::exp(x(1)));
  };

  auto df = [](Eigen::Vector2d x) {
    Eigen::Matrix2d df;
    df << std::exp(x(0)), 0, 0, std::exp(x(1));
    return df;
  };

  // For exp(y): f(y) =  f'(y) for all y and one Newton step satisfies
  // y_{n+1} = y_n - 1
  Eigen::Vector2d x_2ref(-3.0, -1.0);

  // Compute error
  Eigen::Vector2d x_2 = Newton2Steps(f, df, x);
  double err = (x_2 - x_2ref).norm();
  EXPECT_NEAR(0.0, err, 1E-7);
}

}  // namespace MIRK::test