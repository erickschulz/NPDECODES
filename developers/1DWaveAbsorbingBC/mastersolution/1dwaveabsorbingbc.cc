/**
 * @file 1dwaveabsorbingbc.cc
 * @brief NPDE homework "1DWaveAbsorbingBC" code
 * @author Oliver Rietmann
 * @date 08.04.2019
 * @copyright Developed at ETH Zurich
 */

#include "1dwaveabsorbingbc.h"

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
#include <cmath>
#include <utility>
#include <vector>

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

/* SAM_LISTING_BEGIN_1 */
Eigen::MatrixXd waveLeapfrogABC(double c, double T, unsigned int N,
                                unsigned int m) {
  // N is also the number of cells of the mesh
  double h = 1.0 / N;
  // Obtain Galerkin matrices for truncated finite element space
  Eigen::SparseMatrix<double> A = getA_full(N, c, h).block(0, 0, N, N);
  Eigen::SparseMatrix<double> B = getB_full(N).block(0, 0, N, N);
  Eigen::SparseMatrix<double> M = getM_full(N, h).block(0, 0, N, N);
  // Matrix for returning solution
  Eigen::MatrixXd R(m + 1, N + 1);
#if SOLUTION
  Eigen::VectorXd mu = Eigen::VectorXd::Zero(N);   // = mu^(0)
  Eigen::VectorXd nu = Eigen::VectorXd::Zero(N);   // = nu^(-1/2)
  Eigen::VectorXd phi = Eigen::VectorXd::Zero(N);  // = phi(t_0)
  // Universally zero initial conditions make it possible to skip
  // the special initial step usually required for leapfrog.
  double tau = T / m;  // Timestep size
  // The diagonal matrix to be "inverted" in each timestep
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver(1.0 / tau * M + 0.5 * B);
  for (int j = 0; j < m; ++j) {
    R.row(j).head(N) = mu.transpose();
    phi(N - 1) = c / h * g(j * tau);
    nu = solver.solve(-A * mu + (1.0 / tau * M - 0.5 * B) * nu + phi);
    mu = mu + tau * nu;
  }
  R.row(m).head(N) = mu.transpose();

  for (int i = 0; i < m + 1; ++i) {
    R(i, N) = g(i * tau);
  }
#else
  //====================
  // Your code goes here
  //====================
#endif
  return R;
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
std::pair<Eigen::VectorXd, Eigen::VectorXd> computeEnergies(
    const Eigen::MatrixXd &full_solution, double c, double tau) {
  int m = full_solution.rows() - 1;
  int N = full_solution.cols() - 1;
  double h = 1.0 / N;

  Eigen::SparseMatrix<double> A = getA_full(N, c, h);
  Eigen::SparseMatrix<double> M = getM_full(N, h);

  Eigen::VectorXd E_pot(m + 1);
  Eigen::VectorXd E_kin(m);

#if SOLUTION
  Eigen::VectorXd mu0;
  Eigen::VectorXd mu1 = full_solution.row(0).transpose();
  for (int j = 0; j < m; ++j) {
    mu0.swap(mu1);
    E_pot(j) = 0.5 * mu0.dot(A * mu0);
    mu1 = full_solution.row(j + 1).transpose();

    Eigen::VectorXd nu = (mu1 - mu0) / tau;
    E_kin(j) = 0.5 * nu.dot(M * nu);
  }
  E_pot(m) = 0.5 * mu1.dot(A * mu1);
#else
  //====================
  // Your code goes here
  //====================
#endif

  return std::make_pair(E_pot, E_kin);
}
/* SAM_LISTING_END_2 */

}  // namespace WaveAbsorbingBC1D
