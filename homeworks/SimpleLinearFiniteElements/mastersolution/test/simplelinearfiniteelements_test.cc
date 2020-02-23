#include <gtest/gtest.h>

#include <Eigen/SparseLU>

#include "../simplelinearfiniteelements.h"

const double pi = 3.1415926535897;

/**
 * @brief test ElmentMatrix_Mass_LFE implementation
 */
TEST(SimpleLinearFiniteElements, TestElementMatrix_Mass_LFE) {
  // check the produced matrix for a fairly standard triangle
  Eigen::Matrix<double, 2, 3> test;
  test << 0, 1, 0, 0, 0, 1;
  Eigen::Matrix3d M;
  M = SimpleLinearFiniteElements::ElementMatrix_Mass_LFE(test);
  
  Eigen::Matrix3d ref_M;
  ref_M << 0.0833333, 0.0416667, 0.0416667,
  			0.0416667, 0.0833333, 0.0416667,
			0.0416667, 0.0416667, 0.0833333;

  double tol = 1e-8;
  ASSERT_NEAR(ref_M.norm(), M.norm(), tol);
}

/**
 * @brief test L2Error implementation
 */
TEST(SimpleLinearFiniteElements, TestL2Error) {
  // read coarsest mesh
  SimpleLinearFiniteElements::TriaMesh2D square_mesh(CURRENT_SOURCE_DIR "/../../meshes/Square1.txt");
  // exact solution
  auto uExact = [](const Eigen::Vector2d& x) {
    return std::cos(2 * pi * x(0)) * std::cos(2 * pi * x(1));
  };
  // source function
  auto f = [](const Eigen::Vector2d& x) {
    return (8.0 * pi * pi + 1) * std::cos(2 * pi * x(0)) *
           std::cos(2 * pi * x(1));
  };

  // assemble galerkin matrix and load vector
  Eigen::SparseMatrix<double> A = SimpleLinearFiniteElements::GalerkinAssembly(
      square_mesh, SimpleLinearFiniteElements::ElementMatrix_LaplMass_LFE);
  Eigen::VectorXd L = SimpleLinearFiniteElements::assemLoad_LFE(square_mesh, f);

  // solve linear system of equations
  Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> >
      solver;
  solver.analyzePattern(A);
  solver.factorize(A);
  Eigen::VectorXd U = solver.solve(L);

  // compare to expected error
  double error = SimpleLinearFiniteElements::L2Error(square_mesh, U, uExact);
  ASSERT_NEAR(error, 0.232547, 0.1);
}

/**
 * @brief test H1Serror implementation
 */
TEST(SimpleLinearFiniteElements, TestH1Serror) {
  // read coarsest mesh
  SimpleLinearFiniteElements::TriaMesh2D square_mesh(CURRENT_SOURCE_DIR "/../../meshes/Square3.txt");
  // exact gradient
  auto gradUExact = [](const Eigen::Vector2d& x) {
    Eigen::Vector2d gradient;
    gradient << -2 * pi * std::sin(2 * pi * x(0)) * std::cos(2 * pi * x(1)),
        -2 * pi * std::cos(2 * pi * x(0)) * std::sin(2 * pi * x(1));
    return gradient;
  };
  // source function
  auto f = [](const Eigen::Vector2d& x) {
    return (8.0 * pi * pi + 1) * std::cos(2 * pi * x(0)) *
           std::cos(2 * pi * x(1));
  };
  // compute galerkin matrix and load vector
  Eigen::SparseMatrix<double> A = SimpleLinearFiniteElements::GalerkinAssembly(
      square_mesh, SimpleLinearFiniteElements::ElementMatrix_LaplMass_LFE);
  Eigen::VectorXd L = SimpleLinearFiniteElements::assemLoad_LFE(square_mesh, f);

  // solve linear system of equations
  Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> >
      solver;
  solver.analyzePattern(A);
  solver.factorize(A);
  Eigen::VectorXd U = solver.solve(L);

  // compare to expected error
  // high tolerance as the procedure is affected by many rounding errors
  double error =
      SimpleLinearFiniteElements::H1Serror(square_mesh, U, gradUExact);
    ASSERT_NEAR(error, 1.32457, 0.1);
}
