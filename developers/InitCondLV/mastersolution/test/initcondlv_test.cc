/**
 * @file initcondlv_test.cc
 * @brief NPDE homework InitCondLV code
 * @author Oliver Rietmann
 * @date 13.03.2021
 * @copyright Developed at ETH Zurich
 */

#include "../initcondlv.h"

#include <utility>

#include <gtest/gtest.h>

#include <Eigen/Core>

namespace InitCondLV::test {

TEST(InitCondLV, PhiAndW) {

  double u0 = 2.8;
  double v0 = 1.5;
  double T = 2.0;
	std::pair<Eigen::Vector2d, Eigen::Matrix2d> result = PhiAndW(u0, v0, T);

	double tol = 1.0e-8;
  double error;

	Eigen::Vector2d resultPhi = result.first;
  Eigen::Vector2d referencePhi(0.194835295250031, 2.25776610450223);
  error = (referencePhi - resultPhi).lpNorm<Eigen::Infinity>();
  ASSERT_NEAR(0.0, error, tol);
	
	Eigen::Matrix2d resultW = result.second;
  Eigen::Matrix2d referenceW;
  referenceW << -0.158907464246509, 0.0637957462475695, -0.121174713456266, -0.61045517781019;
  error = (referenceW - resultW).lpNorm<Eigen::Infinity>();
  ASSERT_NEAR(0.0, error, tol);
}
	
TEST(InitCondLV, findInitCond) {
	Eigen::Vector2d result = findInitCond();

	Eigen::Vector2d reference(3.1098751029156, 2.08097564048345);

	double tol = 1.0e-8;
  double error = (reference - result).lpNorm<Eigen::Infinity>();
  ASSERT_NEAR(0.0, error, tol);
}

}  // namespace InitCondLV::test
