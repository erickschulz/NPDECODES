/**
 * @file taylorode_test.cc
 * @brief NPDE homework TaylorODE
 * @author ?, Philippe Peter
 * @date 24.03.2021
 * @copyright Developed at ETH Zurich
 */

#include "../taylorode.h"

#include <gtest/gtest.h>

#include <Eigen/Core>
#include <vector>

namespace TaylorODE::test {

TEST(PredPreyModel, f) {
  PredPreyModel model(1.0, 2.0, 4.0, 8.0);
  Eigen::Vector2d y(1.0, 2.0);

  Eigen::Vector2d fy_ex = Eigen::Vector2d(-7, 12);
  Eigen::Vector2d fy_model = model.f(y);

  EXPECT_NEAR(fy_ex(0), fy_model(0), 1E-8);
  EXPECT_NEAR(fy_ex(1), fy_model(1), 1E-8);
}

TEST(PredPreyModel, df) {
  PredPreyModel model(1.0, 2.0, 4.0, 8.0);
  Eigen::Vector2d y(1.0, 2.0);
  Eigen::Vector2d z(1.0, 2.0);

  Eigen::Vector2d fy_ex = Eigen::Vector2d(-15, 28);
  Eigen::Vector2d fy_model = model.df(y, z);

  EXPECT_NEAR(fy_ex(0), fy_model(0), 1E-8);
  EXPECT_NEAR(fy_ex(1), fy_model(1), 1E-8);
}

TEST(PredPreyModel, d2f) {
  PredPreyModel model(1.0, 2.0, 4.0, 8.0);
  Eigen::Vector2d y(1.0, 2.0);
  Eigen::Vector2d z(1.0, 2.0);

  Eigen::Vector2d fy_ex = Eigen::Vector2d(-16, 32);
  Eigen::Vector2d fy_model = model.d2f(y, z);

  EXPECT_NEAR(fy_ex(0), fy_model(0), 1E-8);
  EXPECT_NEAR(fy_ex(1), fy_model(1), 1E-8);
}

TEST(SolvePredPReyTaylor, OneStep) {
  PredPreyModel model(1.0, 2.0, 4.0, 8.0);
  double T = 1.0;
  Eigen::Vector2d y0(1.0, 1.0);
  double M = 1;

  auto res = SolvePredPreyTaylor(model, T, y0, M);
  Eigen::Vector2d y1(14, -43);

  EXPECT_EQ(res.size(), 2);
  EXPECT_NEAR((y0 - res[0]).norm(), 0.0, 1E-8);
  EXPECT_NEAR((y1 - res[1]).norm(), 0.0, 1E-8);
}

}  // namespace TaylorODE::test