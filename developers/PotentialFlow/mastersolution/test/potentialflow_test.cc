/**
 * @file potentialflow_test.cc
 * @brief NPDE homework PotentialFlow code
 * @author R. Hiptmair
 * @date January 2023
 * @copyright Developed at SAM, ETH Zurich
 */

#include "../potentialflow.h"

#include <Eigen/src/SparseLU/SparseLU.h>
#include <gtest/gtest.h>

#include <cstddef>
#include <iostream>

namespace PotentialFlow::test {
TEST(PotentialFlow, printsys) {
  // Constant Dirichlet data
  auto g = [](double /*x*/, double /*y*/) -> double { return 1.0; };
  // Set up linear system
  const unsigned int M = 3;
  Eigen::SparseMatrix<double> A{PotentialFlow::initializeA(M)};
  Eigen::VectorXd rhs{PotentialFlow::initializeRHSVector(g, M)};
  // Output linear system
  Eigen::MatrixXd A_dense = A;
  std::cout << "Galerkin matrix:\n" << A_dense << std::endl;
  std::cout << "rhs vector = " << rhs.transpose() << std::endl;
}

TEST(PotentialFlow, linsol) {
  // Affine linear solutions are reproduced
  auto sol = [](double x, double y) -> double { return 2.0 * x + y + 1; };
  // Set up linear system
  const unsigned int M = 31;
  Eigen::SparseMatrix<double> A{PotentialFlow::initializeA(M)};
  Eigen::VectorXd rhs{PotentialFlow::initializeRHSVector(sol, M)};
  // Solve linear system
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  solver.compute(A);
  if (solver.info() != Eigen::Success) {
    std::cerr << "LU decomposition failed" << std::endl;
  }
  Eigen::VectorXd mu = solver.solve(rhs);
  // Determine maximum norm of the error
  const double h = 1.0 / (M + 1);
  double errmax = 0.0;
  for (unsigned int i = 1; i <= M; ++i) {
    for (unsigned int j = 1; i <= M; ++i) {
      const double sol_val = sol(j * h, i * h);
      const std::size_t idx = (i - 1) * M + j - 1;
      const double mu_val = mu[idx];
      errmax = std::max(errmax, std::abs(sol_val - mu_val));
    }
  }
  EXPECT_NEAR(errmax, 0.0, 1.0E-6);
}
} // namespace PotentialFlow::test
