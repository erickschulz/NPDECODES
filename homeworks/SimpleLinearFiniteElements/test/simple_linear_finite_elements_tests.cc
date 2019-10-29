#include <gtest/gtest.h>

#include "../mysolution/simple_linear_finite_elements.h"

#include <Eigen/SparseLU>

const double pi = 3.1415926535897;

/**
 * @brief test ElmentMatrix_Mass_LFE implementation
 */
TEST(SimpleLinearFiniteElements, TestElementMatrix_Mass_LFE) {
  // check the produced matrix for a fairly standard triangle
  Eigen::Matrix<double, 2, 3> input;
  input << 0, 1, 0, 0, 0, 1;
  Eigen::Matrix3d mat;
  mat = SimpleLinearFiniteElements::ElementMatrix_Mass_LFE(input);
  ASSERT_NEAR(mat(0, 0), 1. / 12., 0.00001);
  ASSERT_NEAR(mat(1, 1), 1. / 12., 0.00001);
  ASSERT_NEAR(mat(2, 2), 1. / 12., 0.00001);
  ASSERT_NEAR(mat(0, 1), 1. / 24., 0.00001);
  ASSERT_NEAR(mat(1, 0), 1. / 24., 0.00001);
  ASSERT_NEAR(mat(2, 0), 1. / 24., 0.00001);
  ASSERT_NEAR(mat(0, 2), 1. / 24., 0.00001);
  ASSERT_NEAR(mat(2, 1), 1. / 24., 0.00001);
  ASSERT_NEAR(mat(1, 2), 1. / 24., 0.00001);

  // check for another, slightly less standard triangle to catch extra errors
  Eigen::Matrix<double, 2, 3> second_input;
  input << 0, 0.9, 0.1, 0.2, .1, 0.3;
  Eigen::Matrix3d second_mat;
  second_mat = SimpleLinearFiniteElements::ElementMatrix_Mass_LFE(input);
  ASSERT_NEAR(second_mat(0, 0), 0.0083333333, 0.00001);
}

/**
 * @brief test L2Error implementation
 */
TEST(SimpleLinearFiniteElements, TestL2Error) {
  // read coarsest mesh
  SimpleLinearFiniteElements::TriaMesh2D square_mesh("Square1.txt");
  // exact solution
  auto uExact = [](double x, double y) {
    return std::cos(2 * pi * x) * std::cos(2 * pi * y);
  };
  // source function
  SimpleLinearFiniteElements::FHandle_t f = [](const Eigen::Vector2d& x) {
    return (8.0 * pi * pi + 1) * std::cos(2 * pi * x(0)) *
           std::cos(2 * pi * x(1));
  };

  // assemble galerkin matrix and load vector
  Eigen::SparseMatrix<double> A = SimpleLinearFiniteElements::GalerkinAssembly(
      square_mesh, SimpleLinearFiniteElements::ElementMatrix_LaplMass_LFE);
  Eigen::VectorXd L = SimpleLinearFiniteElements::assemLoad_LFE(
      square_mesh, SimpleLinearFiniteElements::localLoadLFE, f);

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
  SimpleLinearFiniteElements::TriaMesh2D square_mesh("Square3.txt");
  // exact gradient
  auto gradUExact = [](double x, double y) {
    Eigen::Vector2d gradient;
    gradient << -2 * pi * std::sin(2 * pi * x) * std::cos(2 * pi * y),
        -2 * pi * std::cos(2 * pi * x) * std::sin(2 * pi * y);
    return gradient;
  };
  // source function
  SimpleLinearFiniteElements::FHandle_t f = [](const Eigen::Vector2d& x) {
    return (8.0 * pi * pi + 1) * std::cos(2 * pi * x(0)) *
           std::cos(2 * pi * x(1));
  };
  // compute galerkin matrix and load vector
  Eigen::SparseMatrix<double> A = SimpleLinearFiniteElements::GalerkinAssembly(
      square_mesh, SimpleLinearFiniteElements::ElementMatrix_LaplMass_LFE);
  Eigen::VectorXd L = SimpleLinearFiniteElements::assemLoad_LFE(
      square_mesh, SimpleLinearFiniteElements::localLoadLFE, f);

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
