/**
 * @file nlmatode_test.cc
 * @brief NPDE homework NLMatODE code
 * @author Oliver Rietmann
 * @date 13.03.2021
 * @copyright Developed at ETH Zurich
 */

#include "../nlmatode.h"

#include <gtest/gtest.h>

#include <Eigen/Core>

namespace NLMatODE::test {

double T = 1.0;

Eigen::Matrix3d getY0() {
  Eigen::Matrix3d Y0;
  Y0 << 1.0, 1.0, 0.0, 0.0, 3.0, 2.0, 1.0, 5.0, 2.0;
  return Y0;
}

TEST(NLMatODE, matode) {
  Eigen::Matrix3d result = matode(getY0(), T);

  Eigen::Matrix3d reference;
  reference << 1.11716931299301, 1.06478104796037, -0.317650596136641, 0.858539684847355, 5.24544384944894, 2.58701323916895, 0.121829194605173, 2.5202302578692, 1.09839006062383;

  double tol = 1.0e-8;
  double error = (reference - result).lpNorm<Eigen::Infinity>();
  ASSERT_NEAR(0.0, error, tol);
}

TEST(NLMatODE, checkinvariant) {
	bool result = checkinvariant(getY0(), T); // ode45 does not preserve the norm
  ASSERT_TRUE(!result);
}

}  // namespace NLMatODE::test
