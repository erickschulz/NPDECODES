#include "simple_linear_finite_elements.h"

#include <cmath>
#include <iostream>
#include <tuple>

#include <Eigen/SparseLU>

const double pi = 3.1415926535897;

/**
 * @brief solves system and prints H1-semierror, L2 error, the mesh and a
 * surface plot
 * @param mesh: discretisation of the computational domain
 */
/* SAM_LISTING_BEGIN_1 */
std::tuple<Eigen::VectorXd, double, double> solve(const SimpleLinearFiniteElements::TriaMesh2D& mesh) {

  // define the source function f
  SimpleLinearFiniteElements::FHandle_t f = [](const Eigen::Vector2d& x) {
    return (8.0 * pi * pi + 1) * std::cos(2 * pi * x(0)) *
           std::cos(2 * pi * x(1));
  };

  // the exact solution of the linear variational problem
  auto uExact = [](double x, double y) {
    return std::cos(2 * pi * x) * std::cos(2 * pi * y);
  };

  Eigen::VectorXd U;
  double l2error;
  double h1error;

  /* SOLUTION_BEGIN */
  // the gradient of uExact that can be easily analytically computed
  auto gradUExact = [](double x, double y) {
    Eigen::Vector2d gradient;
    gradient << -2 * pi * std::sin(2 * pi * x) * std::cos(2 * pi * y),
        -2 * pi * std::cos(2 * pi * x) * std::sin(2 * pi * y);
    return gradient;
  };

  // Compute the galerkin matrix A and load vector L
  Eigen::SparseMatrix<double> A = SimpleLinearFiniteElements::GalerkinAssembly(
      mesh, SimpleLinearFiniteElements::ElementMatrix_LaplMass_LFE);
  Eigen::VectorXd L = SimpleLinearFiniteElements::assemLoad_LFE(
      mesh, SimpleLinearFiniteElements::localLoadLFE, f);

  // solve the LSE using the sparse LU solver of Eigen
  Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> >
      solver;
  solver.analyzePattern(A);
  solver.factorize(A);
  U = solver.solve(L);

  // compute both the L2 and the H1s error
  l2error = SimpleLinearFiniteElements::L2Error(mesh, U, uExact);
  h1error = SimpleLinearFiniteElements::H1Serror(mesh, U, gradUExact);
  /* SOLUTION_END */

  return std::make_tuple(U, l2error, h1error);
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
int main() {
  SimpleLinearFiniteElements::TriaMesh2D square_mesh("Square4.txt");
  std::cout << "Mesh loaded " << std::endl;
  std::cout << "Mesh info: " << square_mesh.Coordinates.rows() << " vertices, "
            << square_mesh.Elements.rows() << " elements" << std::endl;

  // print both H1 and L2 errors and plot Mesh
  std::tuple<Eigen::VectorXd, double, double> solution = solve(square_mesh);

  std::cout << "L2-error:  " << std::get<1>(solution) << std::endl;
  std::cout << "H1s-error: " << std::get<2>(solution) << std::endl;

  // plot the mesh and the computed surface
  /* SOLUTION_BEGIN */
  square_mesh.plotMesh("mesh.eps");
  square_mesh.plotSurf("surf.eps", std::get<0>(solution));
  /* SOLUTION_END */
}
/* SAM_LISTING_END_2 */
