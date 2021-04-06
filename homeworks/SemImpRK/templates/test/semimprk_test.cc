/**
 * @file semimprk_test.cc
 * @brief NPDE homework SemImpRK code
 * @author Oliver Rietmann
 * @date 04.04.2021
 * @copyright Developed at ETH Zurich
 */

#include "../semimprk.h"

#include <gtest/gtest.h>

#include <Eigen/Core>

namespace SemImpRK::test {

TEST(SemImpRK, solveRosenbrock) {
  auto f = [](Eigen::Vector3d y) -> Eigen::Vector3d {
    return Eigen::Vector3d(y(0) * y(1), y(1) * y(2), y(2) - y(0));
  };
  auto df = [](Eigen::Vector3d y) {
    Eigen::Matrix3d J;
    J << y(1), y(0), 0.0, 0.0, y(2), y(1), -1.0, 0.0, 1.0;
    return J;
  };
  Eigen::Vector3d y0(1.0, 2.0, 3.0);
  unsigned int N = 10;
  double T = 2.0;

  Eigen::Vector3d yT = SemImpRK::solveRosenbrock(f, df, y0, N, T).back();

  Eigen::Vector3d yT_reference(9.96045385183506, 6.18376651950963e-08,
                               -84.2367784891471);

  double tol = 1.0e-6;
  ASSERT_NEAR(0.0, (yT - yT_reference).lpNorm<Eigen::Infinity>(), tol);
}

TEST(SemImpRK, numexpBurgersGodunov) {
  double convergenceRate = SemImpRK::cvgRosenbrock();
  double convergenceRate_reference = 2.0246142474178366;

  double tol = 1.0e-2;
  ASSERT_NEAR(convergenceRate, convergenceRate_reference, tol);
}

}  // namespace SemImpRK::test
