/**
 * @file engquistoshernumericalflux_test_mastersolution.cc
 * @brief NPDE homework "EngquistOsherNumericalFlux" code
 * @author Oliver Rietmann
 * @date 23.04.2019
 * @copyright Developed at ETH Zurich
 */

#include "../engquistoshernumericalflux.h"

#include <gtest/gtest.h>

#include <Eigen/Core>

namespace EngquistOsherNumericalFlux::test {

TEST(EngquistOsherNumericalFlux, EngquistOsherNumFlux) {
  Eigen::Vector4d v = {-3.0, -2.0, 2.0, 0.5};
  Eigen::Vector4d w = {-1.0, 3.0, 5.0, -1.5};
  Eigen::Vector4d flux_ref = {1.54308063481524, 1.0, 3.76219569108363,
                              2.480035580449628};
  Eigen::Vector4d flux;

  for (int k = 0; k < 4; ++k) {
    flux(k) = EngquistOsherNumFlux(v(k), w(k));
  }

  double tol = 1.0e-8;
  EXPECT_NEAR(flux(0), flux_ref(0), tol);  // v < w < 0
  EXPECT_NEAR(flux(1), flux_ref(1), tol);  // v < 0 < w
  EXPECT_NEAR(flux(2), flux_ref(2), tol);  // 0 < v < w
  EXPECT_NEAR(flux(3), flux_ref(3), tol);  // w < 0 < v
}

TEST(EngquistOsherNumericalFlux, solveCP) {
  Eigen::VectorXd ufinal_ref(6);
  ufinal_ref << 0.0, 0.096997499166518, 0.190962447462037, 0.284836944894129,
      0.215163055105871, 0.109037552537963;

  Eigen::VectorXd uinitial(6);
  uinitial << 0.0, 0.1, 0.2, 0.3, 0.2, 0.1;
  Eigen::VectorXd ufinal = solveCP(0.0, 1.0, uinitial, 0.1);

  double error = (ufinal - ufinal_ref).lpNorm<Eigen::Infinity>();
  double tol = 1.0e-8;
  EXPECT_NEAR(error, 0.0, tol);
}

}  // namespace EngquistOsherNumericalFlux::test
