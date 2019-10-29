/**
 * @ file tests.cc
 * @ brief NPDE homework LinearFE1D code
 * @ author Christian Mitsch
 * @ date 03.03.2019
 * @ copyright Developed at ETH Zurich
 */

#include <gtest/gtest.h>
#include "../mysolution/solve_LinearFE1D.h"

// Test the solver functions
TEST(LinearFE1D, solution_testA) {
  auto gamma = [](double x) { return x; };
  auto f = [](double x) { return x; };
  Vector mesh(9);
  mesh << 0.0, 0.12, 0.2, 0.25, 0.5, 0.7, 0.79, 0.80, 1.0;

  Eigen::VectorXd sol_cor(9);
  sol_cor << 0, 0.0183782, 0.029688, 0.0361259, 0.0592791, 0.0566354, 0.0468281,
      0.0453619, 0;

  Eigen::VectorXd sol = LinearFE1D::solveA(mesh, gamma, f);
  for (int i = 0; i < sol.size(); i++) 
    EXPECT_NEAR(sol(i), sol_cor(i), 1e-5);
}

TEST(LinearFE1D, solution_testB) {
  auto alpha = [](double x) { return x; };
  auto f = [](double x) { return x; };
  Vector mesh(9);
  mesh << 0.0, 0.12, 0.2, 0.25, 0.5, 0.7, 0.79, 0.80, 1.0;

  Eigen::VectorXd sol_cor(9);
  sol_cor << 0.1, 0.412824, 0.48503, 0.514233, 0.576841, 0.570645, 0.556138,
      0.554131, 0.5;

  Eigen::VectorXd sol = LinearFE1D::solveB(mesh, alpha, f, 0.1, 0.5);
  for (int i = 0; i < sol.size(); i++) 
    EXPECT_NEAR(sol(i), sol_cor(i), 1e-5);
}
TEST(LinearFE1D, solution_testC) {
  auto alpha = [](double x) { return x; };
  auto gamma = [](double x) { return x; };
  Vector mesh(9);
  mesh << 0.0, 0.12, 0.2, 0.25, 0.5, 0.7, 0.79, 0.80, 1.0;

  Eigen::VectorXd sol_cor(9);
  sol_cor << 2.42215, 2.30215, 2.23596, 2.19856, 2.04132, 1.96425, 1.94288,
      1.941, 1.92054;

  Eigen::VectorXd sol = LinearFE1D::solveC(mesh, alpha, gamma);
  for (int i = 0; i < sol.size(); i++) 
    EXPECT_NEAR(sol(i), sol_cor(i), 1e-5);
}

/*
//TODO: Implement tests for auxillary functions
TEST(LinearFE1D, test_mat_alpha) {

  Vector mesh(9);
  mesh << 0.0 , 0.12 , 0.2, 0.25, 0.5, 0.7, 0.79, 0.80 ,1.0;
  auto alpha = [](double x) { return x; };

  std::vector<Triplet> alpha_triplets = LinearFE1D::mat_alpha(mesh, alpha);
  SparseMatrix A(9,9);
  A.setFromTriplets(alpha_triplets.begin(), alpha_triplets.end());
  std::cout << A << std::endl;
  EXPECT_EQ(0,42) << "Test for correct solution";
}

TEST(LinearFE1D, test_mat_gamma) {

  Vector mesh(9);
  mesh << 0.0 , 0.12 , 0.2, 0.25, 0.5, 0.7, 0.79, 0.80 ,1.0;
  auto gamma = [](double x) { return x; };

  std::vector<Triplet> gamma_triplets = LinearFE1D::mat_gamma(mesh, gamma);
  //maybe build the sparse matrix for easier testing

  EXPECT_EQ(0,42) << "Test for correct solution";
}

TEST(LinearFE1D, test_rhs_f) {

  Vector mesh(9);
  mesh << 0.0 , 0.12 , 0.2, 0.25, 0.5, 0.7, 0.79, 0.80 ,1.0;
  auto f = [](double x) { return x; };

  Eigen::VectorXd rhs_vector = LinearFE1D::rhs_f(mesh, f);

  EXPECT_EQ(0,42) << "Test for correct solution";
}

TEST(LinearFE1D, test_rhs_constant) {

  Vector mesh(9);
  mesh << 0.0 , 0.12 , 0.2, 0.25, 0.5, 0.7, 0.79, 0.80 ,1.0;
  auto f = [](double x) { return x; };

  Eigen::VectorXd rhs_vector = LinearFE1D::rhs_constant(mesh);

  EXPECT_EQ(0,42) << "Test for correct solution";
}
*/
