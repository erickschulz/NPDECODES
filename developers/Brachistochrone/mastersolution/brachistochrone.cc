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

Eigen::VectorXd coeff_sigma(
    const Eigen::Matrix<double, 2, Eigen::Dynamic> &knots) {
  int M = knots.cols() - 1;
  Eigen::VectorXd sigma(M);

  // Definition quadrature points (2-point Gauss)
  double rh = .5 + .5 / std::sqrt(3.);
  double lh = .5 - .5 / std::sqrt(3.);
  double h = 1. / M;

  // Loop over all segments of the curve
  for (int i = 0; i < M; ++i) {
    // Evaluate the y-component of uh on both quadrature points
    double uh2l = (1 - lh) * knots(1, i) + lh * knots(1, i + 1);
    double uh2r = (1 - rh) * knots(1, i) + rh * knots(1, i + 1);

    // Compute sigma integrated over the segment normalized to unit length
    sigma(i) = 1. / (2. * std::sqrt(-uh2r) * 1. / h *
                     (knots.col(i + 1) - knots.col(i)).lpNorm<2>()) +
               1. / (2. * std::sqrt(-uh2l) * 1. / h *
                     (knots.col(i + 1) - knots.col(i)).lpNorm<2>());
  }

  return sigma;
}

Eigen::VectorXd sourcefn2(
    const Eigen::Matrix<double, 2, Eigen::Dynamic> &knots) {
  int M = knots.cols() - 1;
  Eigen::VectorXd f(M - 1);

  // Definition quadrature points (2-point Gauss)
  double rh = .5 + .5 / std::sqrt(3.);
  double lh = .5 - .5 / std::sqrt(3.);
  double h = 1. / (M);

  // Loop over all segments of the curve
  for (int i = 0; i < M - 1; ++i) {
    // Compute the derivative of uh in the segment 
    // Note that we take the derivative on the right for stability reasons
    Eigen::VectorXd uhp = 1. / h * (knots.col(i + 2) - knots.col(i + 1));

    // Evaluate the y-component of uh on both quadrature points
    double uh2l = (1 - lh) * knots(1, i) + lh * knots(1, i + 1);
    double uh2r = (1 - rh) * knots(1, i) + rh * knots(1, i + 1);

    // Compute the source term
    f(i) = uhp.norm() / (4. * std::sqrt(-uh2r) * uh2r) +
           uhp.norm() / (4. * std::sqrt(-uh2l) * uh2l);
  }

  return f;
}

Eigen::SparseMatrix<double> matR(
    const Eigen::Matrix<double, 2, Eigen::Dynamic> &knots) {
  int M = knots.cols() - 1;
  double h = 1. / (M);

  // Reserve enough data for the matrix
  Eigen::SparseMatrix<double> R(M - 1, M - 1);
  R.reserve(Eigen::RowVectorXi::Constant(M - 1, 3));

  // Compute sigma
  Eigen::VectorXd sigma = coeff_sigma(knots);

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

  // Compress the matrix
  R.makeCompressed();

  return R;
}

Eigen::VectorXd compute_rhs(
    const Eigen::Matrix<double, 2, Eigen::Dynamic> &knots, Eigen::Vector2d b) {
  int M = knots.cols() - 1;
  double h = 1. / M;

  // Compute sigma and source term
  Eigen::VectorXd sigma = coeff_sigma(knots);
  Eigen::VectorXd f = sourcefn2(knots);

  Eigen::VectorXd rhs(2 * (M - 1));

  // Compute first half of rhs
  rhs.setZero();
  rhs(M - 2) = sigma(M - 1) * b(0) / h;

  // Compute second half of rhs
  rhs(M - 1) = sigma(0) * (-0.) / h + h * (f(0));
  for (int i = 1; i < M - 2; ++i) {
    rhs(M - 1 + i) = h * (f(i));
  }
  rhs(2 * M - 3) = sigma(M - 1) * b(1) / h + h * (f(M - 2));

  return rhs;
}

void tabiterr() {
  // Set amount of cells in the mesh
  int M = 128;

  // Initialize b to agree with cycloid
  Eigen::Vector2d b({M_PI, -2.});

  Eigen::MatrixXd knots(2, M + 1);
  for (int i = 0; i < M + 1; ++i) knots(0, i) = (b(0) * i) / M;
  for (int i = 0; i < M + 1; ++i) knots(1, i) = (b(1) * i) / M;

  // Compute and track the brachistochrone function
  std::vector<Eigen::Matrix<double, 2, Eigen::Dynamic>> iterates;
  Eigen::Matrix<double, 2, Eigen::Dynamic> mu = brachistochrone(
      M, b, knots, 1e-12, 1e-12, 1000,
      [&iterates](const Eigen::Matrix<double, 2, Eigen::Dynamic> &mu) -> void {
        iterates.push_back(mu);
      });

  // Compute the L2 norm using 2-point Gauss quadrature
  auto L2norm = [](Eigen::Matrix<double, 2, Eigen::Dynamic> knots) -> double {
    double norm = 0.;

    // Definition quadrature points (2-point Gauss)
    double rh = .5 + .5 / std::sqrt(3.);
    double lh = .5 - .5 / std::sqrt(3.);
    double h = 1. / (knots.cols() - 1.);

    // Loop over all segmenets of the curve
    for (int i = 0; i < knots.cols() - 1; ++i) {
      // Evaluate the y-component of uh on both quadrature points
      Eigen::Vector2d uhl = (1 - lh) * knots.col(i) + lh * knots.col(i + 1);
      Eigen::Vector2d uhr = (1 - rh) * knots.col(i) + rh * knots.col(i + 1);

      // Compute the length associated to this segment of the curve
      norm += h / 2. * (uhl.squaredNorm() + uhr.squaredNorm());
    }
    return std::sqrt(norm);
  };

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

  return;
}

void brachistochrone_cvg() {
  // Set amount of cells in the mesh
  int M = 8;

  // Compute the L2 norm using 2-point Gauss quadrature
  auto L2norm = [](Eigen::Matrix<double, 2, Eigen::Dynamic> knots) -> double {
    double norm = 0.;

    // Definition quadrature points (2-point Gauss)
    double rh = .5 + .5 / std::sqrt(3.);
    double lh = .5 - .5 / std::sqrt(3.);
    double h = 1. / (knots.cols() - 1.);

    // Loop over all segmenets of the curve
    for (int i = 0; i < knots.cols() - 1; ++i) {
      // Evaluate the y-component of uh on both quadrature points
      Eigen::Vector2d uhl = (1 - lh) * knots.col(i) + lh * knots.col(i + 1);
      Eigen::Vector2d uhr = (1 - rh) * knots.col(i) + rh * knots.col(i + 1);

      // Compute the length associated to this segment of the curve
      norm += h / 2. * (uhl.squaredNorm() + uhr.squaredNorm());
    }
    return std::sqrt(norm);
  };

  // Set b to agree with the cycloid
  Eigen::Vector2d b({M_PI, -2.});

  // Initialize parameters for brachistochrone function
  double atol = 1e-6;
  double rtol = 1e-6;
  double itmax = 1000;
  std::vector<double> L2errors;

  // Refine the mesh 8 times and compute the associated L2 errors
  for (int j = 0; j < 8; ++j) {
    // Compute the exact solution (i.e. the cycloid)
    Eigen::VectorXd xi(M + 1);
    Eigen::MatrixXd mu_cycloid(2, M + 1);
    for (int i = 0; i < M + 1; ++i) {
      xi(i) = ((double)i) / M;
      mu_cycloid(0, i) = M_PI * xi(i) - std::sin(M_PI * xi(i));
      mu_cycloid(1, i) = std::cos(M_PI * xi(i)) - 1.;
    }

    // Initialize the knots through a linear approximation
    Eigen::MatrixXd knots(2, M + 1);
    for (int i = 0; i < M + 1; ++i) knots(0, i) = (b(0) * i) / M;
    for (int i = 0; i < M + 1; ++i) knots(1, i) = (b(1) * i) / M;

    // Solve the brachistochrone problem
    Eigen::MatrixXd mu = brachistochrone(M, b, knots, atol, rtol, itmax);

    // Compute and store the L2 error
    L2errors.push_back(L2norm(mu - mu_cycloid));

    // Half the mesh width
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

}  // namespace Brachistochrone
