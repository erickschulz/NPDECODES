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
  EXPECT_NEAR(rhs(2), 2.99022144037, tol);
  EXPECT_NEAR(rhs(3), -5.857963122, tol);
  EXPECT_NEAR(rhs(4), -0.7666687896809623, tol);
  EXPECT_NEAR(rhs(5), -0.87476830362514462, tol);
}

TEST(Brachistochrone, brachistochrone) {
  // We set b to correspond to the cycloid curve
  Eigen::Vector2d b = {M_PI, -2.};

  // We use the following arguments for solving the baristochrone problem
  int M = 80;
  double atol = 1e-10;
  double rtol = 1e-10;
  unsigned itmax = 10000;

  // To make things as easy as possible, we initialize with the exact solution
  Eigen::MatrixXd knots(2, M + 1);
  for (int i = 0; i < M + 1; ++i)
    knots(0, i) = (M_PI * i) / M - std::sin((M_PI * i) / M);
  for (int i = 0; i < M + 1; ++i) knots(1, i) = std::cos((M_PI * i) / M) - 1.;

  // We solve the brachistochrone problem
  Eigen::MatrixXd mu =
      Brachistochrone::brachistochrone(M, b, knots, atol, rtol, itmax);

  // Compute the L2 norm using 2-point Gauss quadrature
  auto L2norm = [](Eigen::Matrix<double, 2, Eigen::Dynamic> knots) -> double {
    double norm = 0.;

    // Definition quadrature points (2-point Gauss)
    double rh = .5 + .5 / std::sqrt(3.);
    double lh = .5 - .5 / std::sqrt(3.);
    double h = 1. / (knots.cols() - 1.);

    // Loop over all segmenets of the curve
    for (int i = 0; i < knots.cols() - 1; ++i) {
      // Evaluate the y-component of uh on both quadrature points
      Eigen::Vector2d uhl = (1 - lh) * knots.col(i) + lh * knots.col(i + 1);
      Eigen::Vector2d uhr = (1 - rh) * knots.col(i) + rh * knots.col(i + 1);

      // Compute the length associated to this segment of the curve
      norm += h / 2. * (uhl.squaredNorm() + uhr.squaredNorm());
    }
    return std::sqrt(norm);
  };

  // We compute the L2 error
  EXPECT_LE(L2norm(mu-knots),1e-1);
}

}  // namespace Brachistochrone::test
