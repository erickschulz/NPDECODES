/**
 * @ file LinearFE1D.cc
 * @ brief NPDE homework LinearFE1D code
 * @ author Christian Mitsch, Amélie Loher
 * @ date 11.11.2019
 * @ copyright Developed at ETH Zurich
 */

// Standard Library
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>
// Eigen
#include <Eigen/SparseCholesky>

namespace LinearFE1D {

typedef Eigen::SparseMatrix<double> SparseMatrix;
typedef Eigen::Triplet<double> Triplet;

// Calculate the matrix entries corresponding to the integral containg alpha
// using composite midpoint rule
/* SAM_LISTING_BEGIN_1 */
template <typename FUNCTOR1>
std::vector<Triplet> computeALap(const Eigen::VectorXd& mesh,
                                 FUNCTOR1&& alpha) {
  unsigned N = mesh.size() - 1;

  std::vector<Triplet> triplets;
  triplets.reserve(3 * N + 1);

  // BEGIN-SOLUTION
  double diag, off_diag;
  double dx_left, dx_right;

  // diagonal entries
  // first and last entry need to be calculated separately
  dx_right = mesh(1) - mesh(0);
  diag = alpha((mesh(1) + mesh(0)) / 2.) / dx_right;
  triplets.push_back(Triplet(0, 0, diag));

  dx_left = mesh(N) - mesh(N - 1);
  diag = alpha((mesh(N) + mesh(N - 1)) / 2.) / dx_left;
  triplets.push_back(Triplet(N, N, diag));

  // other diagonal entries
  for (unsigned i = 1; i < N; ++i) {
    dx_left = mesh(i) - mesh(i - 1);
    dx_right = mesh(i + 1) - mesh(i);
    diag = alpha((mesh(i - 1) + mesh(i)) / 2.) / dx_left +
           alpha((mesh(i) + mesh(i + 1)) / 2.) / dx_right;
    triplets.push_back(Triplet(i, i, diag));
  }

  // off-diagonal entries
  for (unsigned i = 0; i < N; ++i) {
    dx_right = mesh(i + 1) - mesh(i);
    off_diag = -alpha((mesh(i) + mesh(i + 1)) / 2.) / dx_right;
    triplets.push_back(Triplet(i + 1, i, off_diag));
    triplets.push_back(Triplet(i, i + 1, off_diag));
  }
  // END-SOLUTION
  return triplets;
}
/* SAM_LISTING_END_1 */

// Calculate the matrix entries corresponding to the integral containg gamma
// trapezoidal rule
/* SAM_LISTING_BEGIN_2 */
template <typename FUNCTOR1>
std::vector<Triplet> computeMassMatrix(const Eigen::VectorXd& mesh,
                                       FUNCTOR1&& gamma) {
  unsigned N = mesh.size() - 1;

  std::vector<Triplet> triplets;
  triplets.reserve(2 * N + 1);

  // BEGIN-SOLUTION
  double diag, off_diag;
  double dx;

  // diagonal entries
  // first and last entry need to be calculated separately
  dx = mesh(1) - mesh(0);
  diag = gamma(mesh(0)) * 0.5 * dx;
  triplets.push_back(Triplet(0, 0, diag));

  dx = mesh(N) - mesh(N - 1);
  diag = gamma(mesh(N)) * 0.5 * dx;
  triplets.push_back(Triplet(N, N, diag));

  // other diagonal entries
  for (unsigned i = 1; i < N; ++i) {
    dx = mesh(i + 1) - mesh(i - 1);
    diag = gamma(mesh(i)) * 0.5 * dx;
    triplets.push_back(Triplet(i, i, diag));
  }

  // as we are using the trapezoidal rule we do not get any off-diagonal entries

  // END-SOLUTION
  return triplets;
}
/* SAM_LISTING_END_2 */

// Calculate the rhs vector corresponding to the integral containg f*v using
// composite trapezoidal rule
/* SAM_LISTING_BEGIN_3 */
template <typename FUNCTOR1>
Eigen::VectorXd rhs_f(const Eigen::VectorXd& mesh, FUNCTOR1&& f) {
  double dx;
  unsigned N = mesh.size() - 1;
  Eigen::VectorXd b = Eigen::VectorXd::Zero(N + 1);
  // BEGIN-SOLUTION

  // first and last entry need to be calculated separately
  b(0) = f(mesh(0)) * (mesh(1) - mesh(0)) * 0.5;
  b(N) = f(mesh(N)) * (mesh(N) - mesh(N - 1)) * 0.5;

  // other entries
  for (unsigned i = 1; i < N; ++i) {
    dx = mesh(i + 1) - mesh(i - 1);
    b(i) = f(mesh(i)) * 0.5 * dx;
  }

  // END-SOLUTION
  return b;
}

// clang-format on

// Calculate the rhs vector corresponding to the integral containg just v
/* SAM_LISTING_BEGIN_4 */
Eigen::VectorXd rhs_constant(const Eigen::VectorXd& mesh) {
  double dx;
  unsigned N = mesh.size() - 1;
  Eigen::VectorXd b = Eigen::VectorXd::Zero(N + 1);
  // BEGIN-SOLUTION

  // first and last entry need to be calculated separately
  b(0) = (mesh(1) - mesh(0)) / 2.;
  b(N) = (mesh(N) - mesh(N - 1)) / 2.;
  // other entries
  for (unsigned i = 1; i < N; ++i) {
    dx = mesh(i + 1) - mesh(i - 1);
    b(i) = dx;
  }

  // END-SOLUTION
  return b;
}
/* SAM_LISTING_END_4 */

// Build and solve the LSE corresponding to (A)
/* SAM_LISTING_BEGIN_A */
template <typename FUNCTOR1, typename FUNCTOR2>
Eigen::VectorXd solveA(const Eigen::VectorXd& mesh, FUNCTOR1&& gamma,
                       FUNCTOR2&& f) {
  unsigned N = mesh.size() - 1;

  Eigen::VectorXd u(N + 1);
  Eigen::VectorXd b = Eigen::VectorXd::Zero(N + 1);
  // Matrix corresponding to integral with alpha
  Eigen::SparseMatrix<double> A_lap(N + 1, N + 1);
  // Matrix corresponding to integral with gamma
  Eigen::SparseMatrix<double> M(N + 1, N + 1);
  // complete Galerkin Matrix
  Eigen::SparseMatrix<double> A(N + 1, N + 1);

  std::vector<Triplet> triplets_A_Lap;
  std::vector<Triplet> triplets_M;

  /// STEP1: Build the Galerkin matrix A
  /* SOLUTION_BEGIN */

  auto alpha = [](double x) { return 1; };

  triplets_A_Lap = computeALap(mesh, alpha);
  triplets_M = computeMassMatrix(mesh, gamma);

  A_lap.setFromTriplets(triplets_A_Lap.begin(), triplets_A_Lap.end());
  M.setFromTriplets(triplets_M.begin(), triplets_M.end());

  A = A_lap + M;

  /* SOLUTION_END */

  /// STEP2: Build the RHS vector b
  /* SOLUTION_BEGIN */

  b = rhs_f(mesh, f);

  /* SOLUTION_END */

  /// STEP3: Enforce dirichlet boundary conditions
  /* SOLUTION_BEGIN */

  // u(0) = u(N) = 0
  // Thus, we drop first row and column, and last row and column of A and b
  Eigen::SparseMatrix<double> A_reduced(N - 2, N - 2);
  A_reduced = A.block(1, 1, N - 1, N - 1);

  Eigen::VectorXd b_reduced(N - 2);
  b_reduced = b.segment(1, N - 1);

  /* SOLUTION_END */

  /// STEP4: Solve the system Au = b
  /* SOLUTION_BEGIN */
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>, Eigen::Lower> solver;
  solver.compute(A_reduced);

  if (solver.info() != Eigen::Success) {
    throw std::runtime_error("Could not decompose the matrix");
  }

  u(0) = 0;
  u(N) = 0;
  u.segment(1, N - 1) = solver.solve(b_reduced);

  /* SOlUTION_END */

  return u;
}
/* SAM_LISTING_END_A */

// Build an solve the LSE corresponding to (B)
/* SAM_LISTING_BEGIN_B */
template <typename FUNCTOR1, typename FUNCTOR2>
Eigen::VectorXd solveB(const Eigen::VectorXd& mesh, FUNCTOR1&& alpha,
                       FUNCTOR2&& f, double u0, double u1) {
  unsigned N = mesh.size() - 1;
  Eigen::VectorXd u(N + 1);

  double dx_left, dx_right;

  Eigen::VectorXd b = Eigen::VectorXd::Zero(N + 1);
  Eigen::SparseMatrix<double> A(N + 1, N + 1);
  std::vector<Triplet> triplets;

  /// STEP1: Build the Galerkin matrix A
  /* SOLUTION_BEGIN */

  triplets = computeALap(mesh, alpha);
  A.setFromTriplets(triplets.begin(), triplets.end());

  /* SOLUTION_END */

  /// STEP2: Build the RHS vector b
  /* SOLUTION_BEGIN */

  b = rhs_f(mesh, f);

  /* SOLUTION_END */

  /// STEP3: Enforce dirichlet boundary conditions
  /* SOLUTION_BEGIN */

  // Again, we drop first row and column, and last row and column of A and b
  Eigen::SparseMatrix<double> A_reduced(N - 2, N - 2);
  A_reduced = A.block(1, 1, N - 1, N - 1);

  Eigen::VectorXd b_reduced(N - 2);
  b_reduced = b.segment(1, N - 1);

  // Change A and b to enforce non-homogeneous dirichlet boundary contions using
  // the offset function technique
  dx_left = mesh(1) - mesh(0);
  dx_right = mesh(N) - mesh(N - 1);
  b_reduced(0) += u0 * alpha((mesh(0) + mesh(1)) / 2.) / dx_left;
  b_reduced(N - 2) += u1 * alpha((mesh(N - 1) + mesh(N)) / 2.) / dx_right;

  /* SOLUTION_END */

  /// STEP4: Solve the system Au = b
  /* SOLUTION_BEGIN */
  Eigen::SimplicialLDLT<SparseMatrix> solver;
  solver.compute(A_reduced);
  if (solver.info() != Eigen::Success) {
    throw std::runtime_error("Could not decompose the matrix");
  }

  u(0) = u0;
  u(N) = u1;
  u.segment(1, N - 1) = solver.solve(b_reduced);

  /* SOLUTION_END */

  return u;
}
/* SAM_LISTING_END_B */

// Build an solve the LSE corresponding to (C)
/* SAM_LISTING_BEGIN_C */
template <typename FUNCTOR1, typename FUNCTOR2>
Eigen::VectorXd solveC(const Eigen::VectorXd& mesh, FUNCTOR1&& alpha,
                       FUNCTOR2&& gamma) {
  unsigned N = mesh.size() - 1;
  Eigen::VectorXd u(N + 1);
  Eigen::VectorXd b = Eigen::VectorXd::Zero(N + 1);

  Eigen::SparseMatrix<double> A(N + 1, N + 1);
  // Matrix corresponding to integral with alpha
  Eigen::SparseMatrix<double> A_lap(N + 1, N + 1);
  // Matrix corresponding to integral with gamma
  Eigen::SparseMatrix<double> M(N + 1, N + 1);

  std::vector<Triplet> triplets_A_Lap;
  std::vector<Triplet> triplets_M;

  /// STEP1: Build the Galerkin matrix A
  /* SOLUTION_BEGIN */

  triplets_A_Lap = computeALap(mesh, alpha);
  triplets_M = computeMassMatrix(mesh, gamma);
  A_lap.setFromTriplets(triplets_A_Lap.begin(), triplets_A_Lap.end());
  M.setFromTriplets(triplets_M.begin(), triplets_M.end());

  A = A_lap + M;

  /* SOLUTION_END */

  /// STEP2: Build the RHS vector b
  /* SOLUTION_BEGIN */

  b = rhs_constant(mesh);

  // Since we have newmann boundary conditions we don't need to adapt the LSE to
  // enforce boundary conditions
  /* SOLUTION_END */

  /// STEP3: Solve the system Au = b
  /* SOLUTION_BEGIN */
  Eigen::SimplicialLDLT<SparseMatrix> solver;
  solver.compute(A);
  if (solver.info() != Eigen::Success) {
    throw std::runtime_error("Could not decompose the matrix");
  }
  u = solver.solve(b);
  /* SOLUTION_END */

  return u;
}
/* SAM_LISTING_END_C */

}  // namespace LinearFE1D
