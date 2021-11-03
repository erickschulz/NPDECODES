/**
 * @ file master_tests.cc
 * @ brief NPDE homework LinearFE1D code
 * @ author Christian Mitsch, Amélie Loher
 * @ date 11.11.2019
 * @ copyright Developed at ETH Zurich
 */

#include "../linearfe1d.h"

#include <gtest/gtest.h>

#include <Eigen/Core>
#include <Eigen/SparseCore>

// Test the solver functions
TEST(LinearFE1D, solution_testA) {
  auto gamma = [](double x) { return x; };
  auto f = [](double x) { return x; };
  Eigen::VectorXd mesh(9);
  mesh << 0.0, 0.12, 0.2, 0.25, 0.5, 0.7, 0.79, 0.80, 1.0;

  Eigen::VectorXd sol_cor(9);
  sol_cor << 0, 0.0183782, 0.029688, 0.0361259, 0.0592791, 0.0566354, 0.0468281,
      0.0453619, 0;

  Eigen::VectorXd sol = LinearFE1D::solveA(mesh, gamma, f);

  // std::cout << "A" << sol << std::endl;

  for (int i = 0; i < sol.size(); i++) EXPECT_NEAR(sol_cor(i), sol(i), 1e-5);
}

TEST(LinearFE1D, solution_testB) {
  auto alpha = [](double x) { return x; };
  auto f = [](double x) { return x; };
  Eigen::VectorXd mesh(9);
  mesh << 0.0, 0.12, 0.2, 0.25, 0.5, 0.7, 0.79, 0.80, 1.0;

  Eigen::VectorXd sol_cor(9);
  sol_cor << 0.1, 0.412824, 0.48503, 0.514233, 0.576841, 0.570645, 0.556138,
      0.554131, 0.5;

  Eigen::VectorXd sol = LinearFE1D::solveB(mesh, alpha, f, 0.1, 0.5);

  // std::cout << "B" << sol << std::endl;
  for (int i = 0; i < sol.size(); i++) EXPECT_NEAR(sol(i), sol_cor(i), 1e-5);
}

TEST(LinearFE1D, solution_testC) {
  auto alpha = [](double x) { return x; };
  auto gamma = [](double x) { return x; };
  Eigen::VectorXd mesh(9);
  mesh << 0.0, 0.12, 0.2, 0.25, 0.5, 0.7, 0.79, 0.80, 1.0;

  Eigen::VectorXd sol_cor(9);
  sol_cor << 2.42215, 2.30215, 2.23596, 2.19856, 2.04132, 1.96425, 1.94288,
      1.941, 1.92054;

  Eigen::VectorXd sol = LinearFE1D::solveC(mesh, alpha, gamma);

  // std::cout << "C" << sol << std::endl;

  for (int i = 0; i < sol.size(); i++) EXPECT_NEAR(sol(i), sol_cor(i), 1e-5);
}

// TODO: Implement tests for auxillary functions
TEST(LinearFE1D, test_mat_alpha) {
  Eigen::VectorXd mesh(9);
  mesh << 0.0, 0.12, 0.2, 0.25, 0.5, 0.7, 0.79, 0.80, 1.0;
  auto alpha = [](double x) { return x; };

  std::vector<Eigen::Triplet<double>> alpha_triplets =
      LinearFE1D::computeA(mesh, alpha);
  Eigen::SparseMatrix<double> A(9, 9);
  A.setFromTriplets(alpha_triplets.begin(), alpha_triplets.end());
  // std::cout << A << std::endl;

  Eigen::SparseMatrix<double> A_cor(9, 9);

  A_cor.insert(0, 0) = 0.5;
  A_cor.insert(0, 1) = -0.5;
  A_cor.insert(1, 1) = 2.5;
  A_cor.insert(1, 0) = -0.5;
  A_cor.insert(1, 2) = -2.;
  A_cor.insert(2, 2) = 6.5;
  A_cor.insert(2, 1) = -2.;
  A_cor.insert(2, 3) = -4.5;
  A_cor.insert(3, 3) = 6;
  A_cor.insert(3, 2) = -4.5;
  A_cor.insert(3, 4) = -1.5;
  A_cor.insert(4, 4) = 4.5;
  A_cor.insert(4, 3) = -1.5;
  A_cor.insert(4, 5) = -3.;
  A_cor.insert(5, 5) = 11.27778;
  A_cor.insert(5, 4) = -3.;
  A_cor.insert(5, 6) = -8.27778;
  A_cor.insert(6, 6) = 87.77778;
  A_cor.insert(6, 5) = -8.27778;
  A_cor.insert(6, 7) = -79.5;
  A_cor.insert(7, 7) = 84;
  A_cor.insert(7, 6) = -79.5;
  A_cor.insert(7, 8) = -4.5;
  A_cor.insert(8, 8) = 4.5;
  A_cor.insert(8, 7) = -4.5;

  /*A_cor << 0.5, -0.5, 0, 0, 0, 0, 0, 0, 0,
            -0.5, 2.5, -2., 0, 0, 0, 0, 0, 0,
             0, -2., 6.5, -4.5, 0, 0, 0, 0, 0,
                         0, 0, -4.5, 6, -1.5, 0, 0, 0, 0, 0,
             0, 0, 0, -1.5, 4.5, -3., 0, 0, 0,
                         0, 0, 0, 0, -3., 11.27778, -8.27778, 0, 0,
             0, 0, 0, 0, 0, -8.27778, 87.77778, -79.5, 0,
                         0, 0, 0, 0, 0, 0, -79.5, 84, -4.5,
             0, 0, 0, 0, 0, 0, 0, -4.5, 4.5;
  */

  double error = (A_cor - A).norm();
  EXPECT_NEAR(0., error, 1e-5);
}

TEST(LinearFE1D, test_mat_gamma) {
  Eigen::VectorXd mesh(9);
  mesh << 0.0, 0.12, 0.2, 0.25, 0.5, 0.7, 0.79, 0.80, 1.0;
  auto gamma = [](double x) { return x; };

  std::vector<Eigen::Triplet<double>> gamma_triplets =
      LinearFE1D::computeM(mesh, gamma);
  // maybe build the sparse matrix for easier testing
  Eigen::SparseMatrix<double> M(9, 9);
  M.setFromTriplets(gamma_triplets.begin(), gamma_triplets.end());
  // std::cout << M << std::endl;

  Eigen::SparseMatrix<double> M_cor(9, 9);

  M_cor.insert(0, 0) = 0.;
  M_cor.insert(1, 1) = 0.012;
  M_cor.insert(2, 2) = 0.013;
  M_cor.insert(3, 3) = 0.0375;
  M_cor.insert(4, 4) = 0.1125;
  M_cor.insert(5, 5) = 0.1015;
  M_cor.insert(6, 6) = 0.0395;
  M_cor.insert(7, 7) = 0.084;
  M_cor.insert(8, 8) = 0.1;

  double error = (M_cor - M).norm();
  EXPECT_NEAR(0., error, 1e-5);
}

TEST(LinearFE1D, test_rhs_f) {
  Eigen::VectorXd mesh(9);
  mesh << 0.0, 0.12, 0.2, 0.25, 0.5, 0.7, 0.79, 0.80, 1.0;
  auto f = [](double x) { return x; };

  Eigen::VectorXd rhs_vector = LinearFE1D::computeRHS(mesh, f);
  // std::cout << rhs_vector << std::endl;

  Eigen::VectorXd rhs_cor(9);
  rhs_cor << 0., 0.012, 0.013, 0.0375, 0.1125, 0.1015, 0.0395, 0.084, 0.1;

  for (int i = 0; i < rhs_vector.size(); ++i) {
    EXPECT_NEAR(rhs_vector(i), rhs_cor(i), 1e-5);
  }
}
