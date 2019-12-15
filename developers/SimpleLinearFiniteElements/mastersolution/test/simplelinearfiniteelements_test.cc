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
  auto uExact = [](double x, double y) {
    return std::cos(2 * pi * x) * std::cos(2 * pi * y);
  };
  // source function
  std::function<double(const Eigen::Vector2d&)> f = [](const Eigen::Vector2d& x) {
    return (8.0 * pi * pi + 1) * std::cos(2 * pi * x(0)) *
           std::cos(2 * pi * x(1));
  };
  
  // exact solution evaluated at vertices
  Eigen::Vector3d uExact_vec;
  for(size_t t = 0; t < square_mesh.Elements.size(); ++t) {
	const auto& indices = square_mesh.Elements[t];
    for(size_t i = 0; i < 3; ++i) {
	  const auto& v = square_mesh.Vertices[indices(i)];
	  uExact_vec(i) = uExact(v[0],v[1]); 
    }
  }

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
  double error = SimpleLinearFiniteElements::L2Error(square_mesh, U, uExact_vec);
  
  double tol = 1e-4;
  ASSERT_NEAR(error, 0.0611362, tol);
}

/**
 * @brief test H1Serror implementation
 */
TEST(SimpleLinearFiniteElements, TestH1Serror) {
  // read coarsest mesh
  SimpleLinearFiniteElements::TriaMesh2D square_mesh(CURRENT_SOURCE_DIR "/../../meshes/Square3.txt");
  
  // exact gradient
  Eigen::Vector2d gradientExact;
  for(size_t t = 0; t < square_mesh.Elements.size(); ++t) {
    const auto& indices = square_mesh.Elements[t];
    for(size_t i = 0; i < 3; ++i) {
	  const auto& v = square_mesh.Vertices[indices(i)];
	  gradientExact << -2 * pi * std::sin(2 * pi * v[0]) * std::cos(2 * pi * v[1]),
        			   -2 * pi * std::cos(2 * pi * v[0]) * std::sin(2 * pi * v[1]);
     }
  }
  // source function
  std::function<double(const Eigen::Vector2d&)> f = [](const Eigen::Vector2d& x) {
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
      SimpleLinearFiniteElements::H1SError(square_mesh, U, gradientExact);
  
  double tol = 1e-4;
  ASSERT_NEAR(error, 2.5651, tol);
}
