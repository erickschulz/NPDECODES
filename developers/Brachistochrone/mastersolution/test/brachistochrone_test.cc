/**
 * @file XXX_test.cc
 * @brief NPDE homework XXX code
 * @author
 * @date
 * @copyright Developed at SAM, ETH Zurich
 */

#include "../brachistochrone.h"

#include <Eigen/src/Core/util/Constants.h>
#include <gtest/gtest.h>

#include <Eigen/Core>
#include <Eigen/SparseCore>

namespace Brachistochrone::test {

TEST(Brachistochrone, coeff_sigma) {
  // Compute a linear curve
  int M = 4;
  Eigen::Vector2d b(12, -2);
  Eigen::MatrixXd knots(2, M + 1);  // = (Eigen::Matrix<double,2,5>() <<
                                    // 0,-1.,-2.,0.,-.5,-1.).finished();
  for (int i = 0; i < 5; ++i) {
    knots(0, i) = (i * b(0)) / M;
    knots(1, i) = (i * b(1)) / M;
  }

  // Compute sigma
  Eigen::VectorXd sigma = Brachistochrone::coeff_sigma(knots);

  // Test the computed sigma values
  double tol = 1e-4;
  EXPECT_NEAR(sigma(0), 0.1918876, tol);
  EXPECT_NEAR(sigma(1), 0.0962708, tol);
  EXPECT_NEAR(sigma(2), 0.0738926, tol);
  EXPECT_NEAR(sigma(3), 0.0622962, tol);
}

TEST(Brachistochrone, sourcefn2) {
  // Compute a linear curve
  int M = 4;
  Eigen::Vector2d b(1.0, -1.0);
  Eigen::MatrixXd knots(2, M + 1);  
  for (int i = 0; i < 5; ++i) {
    knots(0, i) = (i * b(0)) / M;
    knots(1, i) = (i * b(1)) / M;
  }

  // Compute f
  const Eigen::Matrix<double,Eigen::Dynamic,2> f = Brachistochrone::sourcefn2(knots);

  // Test the computed values for f
  double tol = 1e-2;
  /*
  EXPECT_NEAR(f(0), -100.83227, tol);
  EXPECT_NEAR(f(1), -10.048445, tol);
  EXPECT_NEAR(f(2), -4.4632172, tol);
  */
}

TEST(Brachistochrone, matR) {
  // Compute a linear curve
  int M = 4;
  Eigen::Vector2d b(12, -2);
  Eigen::MatrixXd knots(2, 5);  // = (Eigen::Matrix<double,2,5>() <<
                                // 0,-1.,-2.,0.,-.5,-1.).finished();
  for (int i = 0; i < 5; ++i) {
    knots(0, i) = (i * b(0)) / M;
    knots(1, i) = (i * b(1)) / M;
  }

  // Compute R
  Eigen::SparseMatrix<double> R = Brachistochrone::matR(knots);

  // R should have 7 nonzero values
  EXPECT_EQ(R.nonZeros(), 7);

  // Make R dense to allow evaluation
  Eigen::MatrixXd Rd = R.toDense();

  // Test the computed values for R
  double tol = 1e-4;
  EXPECT_NEAR(Rd(0, 0), 1.15263386784, tol);
  EXPECT_NEAR(Rd(0, 1), -0.385083200, tol);
  EXPECT_NEAR(Rd(1, 0), -0.385083200, tol);
  EXPECT_NEAR(Rd(1, 1), 0.6806539339, tol);
  EXPECT_NEAR(Rd(1, 2), -0.2955707330, tol);
  EXPECT_NEAR(Rd(2, 1), -0.2955707330, tol);
  EXPECT_NEAR(Rd(2, 2), 0.5447558530, tol);
}

TEST(Brachistochrone, compute_rhs) {
  // Compute a linear curve
  int M = 4;
  Eigen::Vector2d b(12, -2);
  Eigen::MatrixXd knots(2, 5);  // = (Eigen::Matrix<double,2,5>() <<
                                // 0,-1.,-2.,0.,-.5,-1.).finished();
  for (int i = 0; i < 5; ++i) {
    knots(0, i) = (i * b(0)) / M;
    knots(1, i) = (i * b(1)) / M;
  }

  // Compute rhs
  Eigen::VectorXd rhs = Brachistochrone::compute_rhs(knots, Eigen::Vector2d(0.0,0.0), b);

  // Test the computed values for rhs
  double tol = 1e-4;
  EXPECT_NEAR(rhs(0), 0., tol);
  EXPECT_NEAR(rhs(1), 0., tol);
  EXPECT_NEAR(rhs(2), 2.9902214403777778, tol);
  EXPECT_NEAR(rhs(3), -8.5620467633785076, tol);
  EXPECT_NEAR(rhs(4), -1.6632919175469014, tol);
  EXPECT_NEAR(rhs(5), -1.3571507626160779, tol);
}



}  // namespace Brachistochrone::test
