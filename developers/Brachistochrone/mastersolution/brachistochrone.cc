/**
 * @file XXX.cc
 * @brief NPDE homework XXX code
 * @author
 * @date
 * @copyright Developed at SAM, ETH Zurich
 */

#include "brachistochrone.h"

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
#include <iomanip>
#include <iostream>

namespace Brachistochrone {

/* SAM_LISTING_BEGIN_0 */
  double L2norm(Eigen::Matrix<double, 2, Eigen::Dynamic> knots) {
    double norm = 0.;
    // Definition quadrature points (2-point Gauss)
    const double rh = .5 + .5 / std::sqrt(3.);
    const double lh = .5 - .5 / std::sqrt(3.);
    double h = 1. / static_cast<double>(knots.cols() - 1);
    // Loop over all segments of the curve
    for (int i = 0; i < knots.cols() - 1; ++i) {
      // Evaluate the y-component of uh on both quadrature points
      const Eigen::Vector2d uhl = (1 - lh) * knots.col(i) + lh * knots.col(i + 1);
      const Eigen::Vector2d uhr = (1 - rh) * knots.col(i) + rh * knots.col(i + 1);
      // Compute the contribution to the L2 norm of that current segment
      norm += h / 2. * (uhl.squaredNorm() + uhr.squaredNorm());
    }
    return std::sqrt(norm);
  }
/* SAM_LISTING_END_0 */

/* SAM_LISTING_BEGIN_1 */
Eigen::VectorXd coeff_sigma(
    const Eigen::Matrix<double, 2, Eigen::Dynamic> &knots) {
  int M = knots.cols() - 1;
  Eigen::VectorXd sigma(M);
  // Definition of locations quadrature points (2-point Gauss rule)
  const double rh = .5 + .5 / std::sqrt(3.);
  const double lh = .5 - .5 / std::sqrt(3.);
  const double h = 1. / M;
  // Loop over all cells of the mesh
  for (int i = 0; i < M; ++i) {
    // Evaluate the y-component of uh on both quadrature points
    const double uh2l = (1 - lh) * knots(1, i) + lh * knots(1, i + 1);
    const double uh2r = (1 - rh) * knots(1, i) + rh * knots(1, i + 1);
    // Compute the derivative
    const Eigen::Vector2d duh = 1. / h * (knots.col(i + 1) - knots.col(i));
    // Norm of the derivative vector
    const double duh_norm = duh.norm();
    // Compute sigma integrated over the segment normalized to unit length
    sigma[i] = 0.5 * (1. / (std::sqrt(-uh2r) * duh_norm) +
                      1. / (std::sqrt(-uh2l) * duh_norm));
  }
  return sigma;
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
Eigen::Matrix<double, Eigen::Dynamic, 2> sourcefn2(
    const Eigen::Matrix<double, 2, Eigen::Dynamic> &knots) {
  const std::size_t M = knots.cols() - 1;
  Eigen::Matrix<double, Eigen::Dynamic, 2> f(M, 2);

  // Definition quadrature points (2-point Gauss)
  const double rh = .5 + .5 / std::sqrt(3.);
  const double lh = .5 - .5 / std::sqrt(3.);
  const double h = 1. / static_cast<double>(M);

  // Loop over all segments of the curve
  for (int i = 0; i < M; ++i) {
    // Compute the derivative of uh in the current cell
    const Eigen::VectorXd uhp = 1. / h * (knots.col(i + 1) - knots.col(i));
    const double uhp_norm = uhp.norm();
    // Evaluate the y-component of uh on both quadrature points
    const double uh2l = (1 - lh) * knots(1, i) + lh * knots(1, i + 1);
    const double uh2r = (1 - rh) * knots(1, i) + rh * knots(1, i + 1);
    // Compute the second components of the function $\Vf$
    f(i, 0) = uhp_norm / (2.0 * std::sqrt(-uh2l) * uh2l);
    f(i, 1) = uhp_norm / (2.0 * std::sqrt(-uh2r) * uh2r);
  }
  return f;
}
/* SAM_LISTING_END_2 */

/* SAM_LISTING_BEGIN_3 */
Eigen::SparseMatrix<double> matR(
    const Eigen::Matrix<double, 2, Eigen::Dynamic> &knots) {
  const Eigen::Index M = knots.cols() - 1;
  double h = 1. / static_cast<double>(M);
  // Reserve enough space for the tridiagonal matrix
  Eigen::SparseMatrix<double> R(M - 1, M - 1);
  R.reserve(Eigen::RowVectorXi::Constant(M - 1, 3));
  // Compute the cell values for the coefficient
  const Eigen::VectorXd sigma = coeff_sigma(knots);
  // Fill in the matrix
  R.insert(0, 0) = (sigma(0) + sigma(1)) / h;
  R.insert(0, 1) = -sigma(1) / h;
  for (int i = 1; i < M + 1 - 3; ++i) {
    R.insert(i, i - 1) = -sigma(i) / h;
    R.insert(i, i) = (sigma(i) + sigma(i + 1)) / h;
    R.insert(i, i + 1) = -sigma(i + 1) / h;
  }
  R.insert(M - 2, M - 3) = -sigma(M - 2) / h;
  R.insert(M - 2, M - 2) = (sigma(M - 2) + sigma(M - 1)) / h;
  // Clean the CRS format
  R.makeCompressed();
  return R;
}
/* SAM_LISTING_END_3 */

/* SAM_LISTING_BEGIN_4 */
Eigen::VectorXd compute_rhs(
    const Eigen::Matrix<double, 2, Eigen::Dynamic> &knots, Eigen::Vector2d b) {
  Eigen::Index M = knots.cols() - 1;
  double h = 1. / static_cast<double>(M);
  // Compute the "coefficients" for all cells
  const Eigen::VectorXd sigma = coeff_sigma(knots);
  const Eigen::Matrix<double, Eigen::Dynamic, 2> f = sourcefn2(knots);
  // Right-hand side vector
  Eigen::VectorXd rhs(2 * (M - 1));
  // Compute head part half of rhs: $\vec{\varphibf}_1$
  rhs.setZero();
  rhs(M - 2) = sigma(M - 1) * b(0) / h;
  // Compute second half of rhs based on the 2-point Gauss rule
  const double w1 = .5 - .5 / std::sqrt(3.);
  const double w2 = .5 - .5 / std::sqrt(3.);
  for (int i = 0; i < M - 1; ++i) {
    rhs(M - 1 + i) =
        0.5 * h *
        (w1 * f(i, 0) + w2 * f(i, 1) + w2 * f(i + 1, 0) + w1 * f(i + 1, 1));
  }
  // Take into account offset function
  rhs(2 * M - 3) += sigma(M - 1) * b(1) / h;
  return rhs;
}
/* SAM_LISTING_END_4 */

/* SAM_LISTING_BEGIN_5 */
void tabiterr(unsigned int M, double rtol, double atol, unsigned int maxit) {
  // Initialize b to agree with endpoint of cycloid curve
  Eigen::Vector2d b({M_PI, -2.});
  // Initial guess: straight line
  Eigen::MatrixXd knots(2, M + 1);
  for (int i = 0; i < M + 1; ++i) knots(0, i) = (b(0) * i) / M;
  for (int i = 0; i < M + 1; ++i) knots(1, i) = (b(1) * i) / M;
  // Carry out the fixed-point iteration
  std::vector<Eigen::Matrix<double, 2, Eigen::Dynamic>> iterates;
  Eigen::Matrix<double, 2, Eigen::Dynamic> mu = brachistochrone(
      M, b, knots, atol, rtol, maxit,
      [&iterates](const Eigen::Matrix<double, 2, Eigen::Dynamic> &mu) -> void {
        iterates.push_back(mu);
      });
  // Print table with results
  std::cout << "| k   | ||u_h^(k)-u_h^(" << iterates.size() - 1
            << ")||_{L2} |\n";
  std::cout << "-----------------------------------|\n";
  for (int i = 0; i < iterates.size() - 1; ++i) {
    std::cout << "| " << std::setw(3) << i << " | " << std::setw(26)
              << std::setprecision(5) << std::scientific
              << L2norm(iterates.at(i) - iterates.at(iterates.size() - 1))
              << " |\n";
  }
}
/* SAM_LISTING_END_5 */

/* SAM_LISTING_BEGIN_6 */
void brachistochrone_cvg() {
  // Set minimal number of cells in the mesh
  int M = 8;
  // Set b to the endpoint of the cycloid curve
  Eigen::Vector2d b({M_PI, -2.});
  // Initialize parameters for brachistochrone function
  double atol = 1e-6;
  double rtol = 1e-8;
  double itmax = 10000;
  std::vector<double> L2errors;

  // Refine the mesh 8 times and compute the associated L2 errors
  for (int j = 0; j < 8; ++j) {
    // Compute the exact solution (i.e. the cycloid)
    Eigen::VectorXd xi(M + 1);
    Eigen::MatrixXd mu_cycloid(2, M + 1);
    for (int i = 0; i <= M; ++i) {
      xi(i) = static_cast<double>(i) / M;
      mu_cycloid(0, i) = M_PI * xi(i) - std::sin(M_PI * xi(i));
      mu_cycloid(1, i) = std::cos(M_PI * xi(i)) - 1.;
    }
    // Initialize the knots through a linear approximation
    Eigen::MatrixXd knots(2, M + 1);
    for (int i = 0; i <= M; ++i) knots(0, i) = (b(0) * i) / M;
    for (int i = 0; i <= M; ++i) knots(1, i) = (b(1) * i) / M;
    // Solve the brachistochrone problem
    Eigen::MatrixXd mu = brachistochrone(M, b, knots, atol, rtol, itmax);
    // Compute and store the L2 error
    L2errors.push_back(L2norm(mu - mu_cycloid));
    // Half the mesh width/double the number of cells
    M *= 2;
  }
  // Print a table with results.
  std::cout << "| k   | ||u_h-I_M u_h||_{L2} |\n";
  std::cout << "-----------------------------|\n";
  for (int i = 0; i < L2errors.size() - 1; ++i) {
    std::cout << "| " << std::setw(3) << 8 * std::pow(2, i) << " | "
              << std::setw(20) << std::setprecision(5) << std::scientific
              << L2errors.at(i) << " |\n";
  }
}
/* SAM_LISTING_END_6 */

}  // namespace Brachistochrone
