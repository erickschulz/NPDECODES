/**
 * @file XXX.h
 * @brief NPDE homework XXX code
 * @author
 * @date
 * @copyright Developed at SAM, ETH Zurich
 */

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
#include <iostream>

// #include <lf/assemble/assemble.h>
// #include <lf/fe/fe.h>
// #include <lf/mesh/utils/utils.h>
// #include <lf/uscalfe/uscalfe.h>

namespace Brachistochrone {

Eigen::VectorXd coeff_sigma(
    const Eigen::Matrix<double, 2, Eigen::Dynamic> &knots);

Eigen::Matrix<double, Eigen::Dynamic, 2> sourcefn2(
    const Eigen::Matrix<double, 2, Eigen::Dynamic> &knots);

Eigen::SparseMatrix<double> matR(
    const Eigen::Matrix<double, 2, Eigen::Dynamic> &knots);

Eigen::VectorXd compute_rhs(
    const Eigen::Matrix<double, 2, Eigen::Dynamic> &knots, Eigen::Vector2d b);

double L2norm(Eigen::Matrix<double, 2, Eigen::Dynamic> knots);

/* SAM_LISTING_BEGIN_1 */
template <typename RECORDER = std::function<
              void(const Eigen::Matrix<double, 2, Eigen::Dynamic> &)>>
Eigen::Matrix<double, 2, Eigen::Dynamic> brachistochrone(
    unsigned int M, Eigen::Vector2d b,
    const Eigen::Matrix<double, 2, Eigen::Dynamic> &knots0, double atol,
    double rtol, unsigned int itmax,
    RECORDER &&rec = [](const Eigen::Matrix<double, 2, Eigen::Dynamic> &)
        -> void { return; }) {
  // Initialize knots
  Eigen::MatrixXd knots = knots0;
  Eigen::MatrixXd knots_prev = knots;

  // Computes the length of the discrete curve
  auto lencomp = [&knots]() {
    double out = 0.;

    // Definition quadrature points (2-point Gauss)
    double rh = .5 + .5 / std::sqrt(3.);
    double lh = .5 - .5 / std::sqrt(3.);
    double h = 1. / (knots.cols() - 1);

    // Loop over all segments of the curve
    for (int i = 0; i < knots.cols() - 1; ++i) {
      // Evaluate the y-component of uh on both quadrature points
      double uh2l = (1 - lh) * knots(1, i) + lh * knots(1, i + 1);
      double uh2r = (1 - rh) * knots(1, i) + rh * knots(1, i + 1);

      // Compute the length associated to this segment of the curve
      out += h * (knots.col(i) - knots.col(i + 1)).norm() * .5 *
             (1. / std::sqrt(-uh2l) + 1. / std::sqrt(-uh2r));
    }
    return out;
  };
  // Record the initial conditions of the knots
  rec(knots);
  // Perform the fixed point-iteration
  unsigned int q = 0;
  do {
    // Store the previous knots data for comparison later
    knots_prev = knots;

    // Compute R-matrix and RHS
    Eigen::SparseMatrix<double> R = matR(knots);
    Eigen::VectorXd rhs = compute_rhs(knots, b);
    // Initialize the LU solver
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.compute(R);
    // Solve the first system
    Eigen::VectorXd uh = solver.solve(rhs.segment(0, M - 1));
    for (int i = 0; i < M - 1; ++i) knots(0, i + 1) = uh(i);
    // Solve the second system
    Eigen::VectorXd uh2 = solver.solve(rhs.segment(M - 1, M - 1));
    for (int i = 0; i < M - 1; ++i) knots(1, i + 1) = uh2(i);
    // Update and record
    rec(knots);
    // Abort when max iterations is reached or error is small enough.
  } while ((++q < itmax) &&
           (Brachistochrone::L2norm(knots - knots_prev) >
            std::min<double>(atol, rtol * Brachistochrone::L2norm(knots))));

  // Output a warning if fixed-point iteration was aborted because max
  // iterations was reached.
  if (q == itmax)
    std::cout << "Warning: Max iterations reached. Truncated with L2 error: "
              << L2norm(knots - knots_prev) << "\n";
  return knots;
};

/* SAM_LISTING_END_1 */

void tabiterr(unsigned int M = 128, double rtol = 1.0E-4, double atol = 1.0E-6,
              unsigned int maxit = 10000);

void brachistochrone_cvg(void);

}  // namespace Brachistochrone
