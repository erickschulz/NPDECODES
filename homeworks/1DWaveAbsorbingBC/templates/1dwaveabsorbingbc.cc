/**
 * @file 1dwaveabsorbingbc.cc
 * @brief NPDE homework "1DWaveAbsorbingBC" code
 * @author Oliver Rietmann
 * @date 08.04.2019
 * @copyright Developed at ETH Zurich
 */

#include "1dwaveabsorbingbc.h"

#include <cmath>
#include <utility>
#include <vector>

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>

namespace WaveAbsorbingBC1D {

constexpr double PI = 3.14159265358979323846;

// Boundary data on the right-hand side (x = 1)
double g(double t) { return 0 <= t && t <= PI ? std::sin(t) : 0.0; }

/**
 * @brief Get the full (--> including both boundary points) Galerkin matrix A
 * @param N number of spacial nodes, including x=0, but excluding x=1
 * @param c speed of propagation
 * @param h (spacial) meshwidth
 * @return Full Galerkin matrix A of shape (N+2) times (N+2)
 */
Eigen::SparseMatrix<double> getA_full(unsigned int N, double c, double h) {
  std::vector<Eigen::Triplet<double>> triplets;
  triplets.reserve(3 * (N + 1) - 2);

  double scale = c / h;

  // store first row separately
  triplets.push_back(Eigen::Triplet<double>(0, 0, scale));
  triplets.push_back(Eigen::Triplet<double>(0, 1, -scale));

  // loop over interior rows
  for (unsigned i = 1; i < N; ++i) {
    triplets.push_back(Eigen::Triplet<double>(i, i - 1, -scale));
    triplets.push_back(Eigen::Triplet<double>(i, i, 2.0 * scale));
    triplets.push_back(Eigen::Triplet<double>(i, i + 1, -scale));
  }

  // store last row separately
  triplets.push_back(Eigen::Triplet<double>(N, N - 1, -scale));
  triplets.push_back(Eigen::Triplet<double>(N, N, scale));

  Eigen::SparseMatrix<double> A(N + 1, N + 1);
  A.setFromTriplets(triplets.begin(), triplets.end());
  return A;
}

/**
 * @brief Get the full (--> including both boundary points) Galerkin matrix B
 * @param N number of spacial nodes, including x=0, but excluding x=1
 * @return Full Galerkin matrix B of shape (N+2) times (N+2)
 */
Eigen::SparseMatrix<double> getB_full(unsigned int N) {
  Eigen::SparseMatrix<double> B(N + 1, N + 1);
  B.coeffRef(0, 0) = 1.0;
  return B;
}

/**
 * @brief Get the full (--> including both boundary points) Galerkin matrix M
 * @param N number of spacial nodes, including x=0, but excluding x=1
 * @param h (spacial) meshwidth
 * @return Full Galerkin matrix M of shape (N+2) times (N+2)
 */
Eigen::SparseMatrix<double> getM_full(unsigned int N, double h) {
  Eigen::SparseMatrix<double> M(N + 1, N + 1);
  M.setIdentity();
  M *= h;
  M.coeffRef(0, 0) = h / 2;
  M.coeffRef(N, N) = h / 2;
  return M;
}

Eigen::MatrixXd waveLeapfrogABC(double c, double T, unsigned int N,
                                unsigned int m) {
  double h = 1.0 / N;

  Eigen::SparseMatrix<double> A = getA_full(N, c, h).block(0, 0, N, N);
  Eigen::SparseMatrix<double> B = getB_full(N).block(0, 0, N, N);
  Eigen::SparseMatrix<double> M = getM_full(N, h).block(0, 0, N, N);

  Eigen::MatrixXd R(m + 1, N + 1);

  /* SOLUTION_BEGIN */
  /* SOLUTION_END */

  return R;
}

std::pair<Eigen::VectorXd, Eigen::VectorXd> computeEnergies(
    const Eigen::MatrixXd &full_solution, double c, double tau) {
  int m = full_solution.rows() - 1;
  int N = full_solution.cols() - 1;
  double h = 1.0 / N;

  Eigen::SparseMatrix<double> A = getA_full(N, c, h);
  Eigen::SparseMatrix<double> M = getM_full(N, h);

  Eigen::VectorXd E_pot(m + 1);
  Eigen::VectorXd E_kin(m);

  /* SOLUTION_BEGIN */
  /* SOLUTION_END */

  return std::make_pair(E_pot, E_kin);
}


}  // namespace WaveAbsorbingBC1D
