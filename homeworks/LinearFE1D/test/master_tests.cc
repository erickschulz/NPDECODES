/**
 * @ file master_tests.cc
 * @ brief NPDE homework LinearFE1D code
 * @ author Christian Mitsch
 * @ date 01.03.2019
 * @ copyright Developed at ETH Zurich
 */

#include <gtest/gtest.h>
#include "../mastersolution/linearfe1d.h"

// Test the solver functions
TEST(LinearFE1D, solution_testA) {
  auto gamma = [](double x) { return x; };
  auto f = [](double x) { return x; };
  Vector mesh(9);
  mesh << 0.0, 0.12, 0.2, 0.25, 0.5, 0.7, 0.79, 0.80, 1.0;

  Vector sol_cor(9);
  sol_cor << 0, 0.0183782, 0.029688, 0.0361259, 0.0592791, 0.0566354, 0.0468281,
      0.0453619, 0;

  Vector sol = LinearFE1D::solveA(mesh, gamma, f);
  for (int i = 0; i < sol.size(); i++) 
    EXPECT_NEAR(sol(i), sol_cor(i), 1e-5);
}

TEST(LinearFE1D, solution_testB) {
  auto alpha = [](double x) { return x; };
  auto f = [](double x) { return x; };
  Vector mesh(9);
  mesh << 0.0, 0.12, 0.2, 0.25, 0.5, 0.7, 0.79, 0.80, 1.0;

  Vector sol_cor(9);
  sol_cor << 0.1, 0.412824, 0.48503, 0.514233, 0.576841, 0.570645, 0.556138,
      0.554131, 0.5;

  Vector sol = LinearFE1D::solveB(mesh, alpha, f, 0.1, 0.5);
  for (int i = 0; i < sol.size(); i++) 
    EXPECT_NEAR(sol(i), sol_cor(i), 1e-5);
}
TEST(LinearFE1D, solution_testC) {
  auto alpha = [](double x) { return x; };
  auto gamma = [](double x) { return x; };
  Vector mesh(9);
  mesh << 0.0, 0.12, 0.2, 0.25, 0.5, 0.7, 0.79, 0.80, 1.0;

  Vector sol_cor(9);
  sol_cor << 2.42215, 2.30215, 2.23596, 2.19856, 2.04132, 1.96425, 1.94288,
      1.941, 1.92054;

  Vector sol = LinearFE1D::solveC(mesh, alpha, gamma);
  for (int i = 0; i < sol.size(); i++) 
    EXPECT_NEAR(sol(i), sol_cor(i), 1e-5);
}


//TODO: Implement tests for auxillary functions
TEST(LinearFE1D, test_mat_alpha) {

  Vector mesh(9);
  mesh << 0.0 , 0.12 , 0.2, 0.25, 0.5, 0.7, 0.79, 0.80 ,1.0;
  auto alpha = [](double x) { return x; };

  std::vector<Triplet> alpha_triplets = LinearFE1D::mat_alpha(mesh, alpha);
  SparseMatrix A(9,9);
  A.setFromTriplets(alpha_triplets.begin(), alpha_triplets.end());
  //std::cout << A << std::endl;
  
  SparseMatrix A_cor(9,9); 

  A_cor.insert(0,0) = 0.5; A_cor.insert(0,1) = -0.5;
  A_cor.insert(1,1) = 2.5; A_cor.insert(1,0) = -0.5; A_cor.insert(1,2) = -2.;
  A_cor.insert(2,2) = 6.5; A_cor.insert(2,1) = -2.; A_cor.insert(2,3) = -4.5;
  A_cor.insert(3,3) = 6; A_cor.insert(3,2) = -4.5; A_cor.insert(3,4) = -1.5;
  A_cor.insert(4,4) = 4.5; A_cor.insert(4,3) = -1.5; A_cor.insert(4,5) = -3.;
  A_cor.insert(5,5) = 11.27778; A_cor.insert(5,4) = -3.; A_cor.insert(5,6) = -8.27778;
  A_cor.insert(6,6) = 87.77778; A_cor.insert(6,5) = -8.27778; A_cor.insert(6,7) = -79.5;
  A_cor.insert(7,7) = 84; A_cor.insert(7,6) = -79.5; A_cor.insert(7,8) = -4.5;
  A_cor.insert(8,8) = 4.5; A_cor.insert(8,7) = -4.5;

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

  Vector mesh(9);
  mesh << 0.0 , 0.12 , 0.2, 0.25, 0.5, 0.7, 0.79, 0.80 ,1.0;
  auto gamma = [](double x) { return x; };

  std::vector<Triplet> gamma_triplets = LinearFE1D::mat_gamma(mesh, gamma);
  //maybe build the sparse matrix for easier testing
  SparseMatrix M(9,9);
  M.setFromTriplets(gamma_triplets.begin(), gamma_triplets.end());
  //std::cout << M << std::endl;
  
  SparseMatrix M_cor(9,9);

  M_cor.insert(0,0) = 0.;
  M_cor.insert(1,1) = 0.012;
  M_cor.insert(2,2) = 0.013;
  M_cor.insert(3,3) = 0.0375;
  M_cor.insert(4,4) = 0.1125;
  M_cor.insert(5,5) = 0.1015;
  M_cor.insert(6,6) = 0.0395;
  M_cor.insert(7,7) = 0.084;
  M_cor.insert(8,8) = 0.1;

  double error = (M_cor - M).norm();
  EXPECT_NEAR(0., error, 1e-5);
}

TEST(LinearFE1D, test_rhs_f) {

  Vector mesh(9);
  mesh << 0.0 , 0.12 , 0.2, 0.25, 0.5, 0.7, 0.79, 0.80 ,1.0;
  auto f = [](double x) { return x; };

  Vector rhs_vector = LinearFE1D::rhs_f(mesh, f);
  //std::cout << rhs_vector << std::endl;
  
  Vector rhs_cor(9);
  rhs_cor << 0., 0.012, 0.013, 0.0375, 0.1125, 0.1015, 0.0395, 0.084, 0.1;

  for(int i = 0; i < rhs_vector.size(); ++i) {
	  EXPECT_NEAR(rhs_vector(i), rhs_cor(i), 1e-5);
  }

}

TEST(LinearFE1D, test_rhs_constant) {

  Vector mesh(9);
  mesh << 0.0 , 0.12 , 0.2, 0.25, 0.5, 0.7, 0.79, 0.80 ,1.0;
  auto f = [](double x) { return x; };

  Vector rhs_vector = LinearFE1D::rhs_constant(mesh);
  //std::cout << rhs_vector << std::endl;

  Vector rhs_cor(9);
  rhs_cor << 0.06, 0.2, 0.13, 0.3, 0.45, 0.29, 0.1, 0.21, 0.1;

  for(int i = 0; i < rhs_vector.size(); ++i) {
      EXPECT_NEAR(rhs_vector(i), rhs_cor(i), 1e-5);
  }

}

