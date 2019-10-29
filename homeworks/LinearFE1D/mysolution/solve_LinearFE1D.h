/**
 * @ file LinearFE1D.cc
 * @ brief NPDE homework LinearFE1D code
 * @ author Christian Mitsch
 * @ date 01.03.2019
 * @ copyright Developed at ETH Zurich
 */

#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <cmath>
#include <iostream>
#include <vector>

typedef Eigen::SparseMatrix<double> SparseMatrix;
typedef Eigen::Triplet<double> Triplet;
typedef Eigen::VectorXd Vector;

namespace LinearFE1D {

// Calculate the matrix entries corresponding to the integral containg alpha
template <typename FUNCTOR1>
std::vector<Triplet> mat_alpha(const Eigen::VectorXd &mesh, FUNCTOR1 alpha) {
  std::vector<Triplet> triplets;
  unsigned M = mesh.size() - 1;
  triplets.reserve(3 * M + 1);

  // BEGIN_SOLUTION
  // TODO Your implementation goes here!
  // END_SOLUTION
  return triplets;
}

// Calculate the matrix entries corresponding to the integral containg gamma
template <typename FUNCTOR1>
std::vector<Triplet> mat_gamma(const Eigen::VectorXd &mesh, FUNCTOR1 gamma) {
  std::vector<Triplet> triplets;
  unsigned M = mesh.size() - 1;
  triplets.reserve(M + 1);

  // BEGIN_SOLUTION
  // TODO Your implementation goes here!
  // END_SOLUTION
  return triplets;
}

// Calculate the rhs vector corresponding to the integral containg f*v using
// composite trapezoidal rule
template <typename FUNCTOR1>
Eigen::VectorXd rhs_f(const Eigen::VectorXd &mesh, FUNCTOR1 f) {
  double dx_left, dx_right;
  unsigned M = mesh.size() - 1;
  Vector b = Vector::Zero(M + 1);
  // BEGIN_SOLUTION
  // TODO Your implementation goes here!
  // END_SOLUTION
  return b;
}

// Calculate the rhs vector corresponding to the integral containg just v
Eigen::VectorXd rhs_constant(const Eigen::VectorXd &mesh) {
  double dx_left, dx_right;
  unsigned M = mesh.size() - 1;
  Vector b = Vector::Zero(M + 1);
  // BEGIN_SOLUTION
  // TODO Your implementation goes here!
  // END_SOLUTION
  return b;
}

// Build an solve the LSE corresponding to (A)
template <typename FUNCTOR1, typename FUNCTOR2>
Eigen::VectorXd solveA(const Eigen::VectorXd &mesh, FUNCTOR1 gamma,
                       FUNCTOR2 f) {
  unsigned M = mesh.size() - 1;
  Vector u(M + 1);

  Vector b = Vector::Zero(M + 1);
  SparseMatrix A(M + 1, M + 1);
  std::vector<Triplet> triplets;
  triplets.reserve(2 * M + 1);

  /// STEP1: Build the Galerkin matrix A
  // BEGIN_SOLUTION
  // TODO Your implementation goes here!
  // END_SOLUTION

  /// STEP2: Build the RHS vector b
  // BEGIN_SOLUTION
  // TODO Your implementation goes here!
  // END_SOLUTION

  /// STEP3: Enforce dirichlet boundary conditions
  // BEGIN_SOLUTION
  // TODO Your implementation goes here!
  // END_SOLUTION

  /// STEP4: Solve the system Au = b
  // BEGIN_SOLUTION
  // TODO Your implementation goes here!
  // END_SOLUTION

  return u;
}

// Build an solve the LSE corresponding to (B)
template <typename FUNCTOR1, typename FUNCTOR2>
Eigen::VectorXd solveB(const Eigen::VectorXd &mesh, FUNCTOR1 alpha, FUNCTOR2 f,
                       double u0, double u1) {
  unsigned M = mesh.size() - 1;
  Vector u(M + 1);

  Vector b = Vector::Zero(M + 1);
  SparseMatrix A(M + 1, M + 1);
  std::vector<Triplet> triplets;
  triplets.reserve(2 * M + 1);

  /// STEP1: Build the Galerkin matrix A
  // BEGIN_SOLUTION
  // TODO Your implementation goes here!
  // END_SOLUTION

  /// STEP2: Build the RHS vector b
  // BEGIN_SOLUTION
  // TODO Your implementation goes here!
  // END_SOLUTION

  /// STEP3: Enforce dirichlet boundary conditions
  // BEGIN_SOLUTION
  // TODO Your implementation goes here!
  // END_SOLUTION

  /// STEP4: Solve the system Au = b
  // BEGIN_SOLUTION
  // TODO Your implementation goes here!
  // END_SOLUTION

  return u;
}

// Build an solve the LSE corresponding to (C)
template <typename FUNCTOR1, typename FUNCTOR2>
Eigen::VectorXd solveC(const Eigen::VectorXd &mesh, FUNCTOR1 alpha,
                       FUNCTOR2 gamma) {
  unsigned M = mesh.size() - 1;
  Vector u(M + 1);

  Vector b = Vector::Zero(M + 1);
  SparseMatrix A(M + 1, M + 1);
  std::vector<Triplet> triplets;
  triplets.reserve(2 * M + 1);

  /// STEP1: Build the Galerkin matrix A
  // BEGIN_SOLUTION
  // TODO Your implementation goes here!
  // END_SOLUTION

  /// STEP2: Build the RHS vector b
  // BEGIN_SOLUTION
  // TODO Your implementation goes here!
  // END_SOLUTION

  /// STEP3: Solve the system Au = b
  // BEGIN_SOLUTION
  // TODO Your implementation goes here!
  // END_SOLUTION

  return u;
}
}  // namespace LinearFE1D
