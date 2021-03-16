/**
 * @file systemode_main.cc
 * @brief NPDE homework SystemODE
 * @copyright Developed at ETH Zurich
 */

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>

#include <cmath>
#include <iostream>
#include <iomanip>

#include "systemode.h"

/* SAM_LISTING_BEGIN_0 */
int main() {
  // PARAMETERS
  double T = 1;
  int n = 5;

  // INITIAL VALUE
  Eigen::VectorXd y0(2 * n);
  for (int i = 0; i < n; ++i) {
    y0(i) = (i + 1.) / n;
    y0(i + n) = -1;
  }

  // SETUP
  double conv_rate = 0;
  std::cout << std::setw(8) << "M" << std::setw(20) << "Error" << std::endl;

  // IMPLEMENTATION OF RIGHT-HAND SIDE f
  // Build tridiagonal C matrix
  Eigen::SparseMatrix<double> C(n, n);
  C.reserve(Eigen::VectorXi::Constant(n, 3));
  C.insert(0, 0) = 2;
  for (int i = 1; i < n; ++i) {
    C.insert(i, i) = 2;
    C.insert(i, i - 1) = -1;
    C.insert(i - 1, i) = -1;
  }
  C.makeCompressed();
  // Compute the right-hand side f
  // The system of ODEs is y' = f(y) with y = [u;v],
  // and f(y) = f([u;v]) = [v;C^{-1}r(u)]
  auto f = [n, C](Eigen::VectorXd y) {
    Eigen::VectorXd fy(2 * n);
    fy.head(n) = y.tail(n);
    Eigen::VectorXd r(n);
    r(0) = y(0) * (y(1) + y(0));
    r(n - 1) = y(n - 1) * (y(n - 1) + y(n - 2));
    for (int i = 1; i < n - 1; ++i) {
      r(i) = y(i) * (y(i - 1) + y(i + 1));
    }
    Eigen::SparseLU<Eigen::SparseMatrix<double>> Csolver;
    Csolver.compute(C);
    fy.tail(n) = Csolver.solve(r);
    return fy;
  };

  // COMPUTE AN "EXACT" SOLUTION
  // Use N=2^12 steps to calculate an approximate "exact" solution.
  int N_exact = std::pow(2, 12); // number of steps
  double h = T / N_exact;  // step size
  Eigen::VectorXd yT_exact = y0; // initial value
  Eigen::VectorXd y_next;
  for (int step = 0; step < N_exact; step++) {
    y_next = SystemODE::rk4step(f, h, yT_exact);
    yT_exact = y_next;
  }

  // CONVERGENCE ANALYSIS
  // Calculate solution using N=2,...,2^kmax steps.
  int kmax = 10;
  Eigen::VectorXd Error(kmax);
  for (int k = 0; k < kmax; k++) {
    int M = std::pow(2, k + 1); // number of steps
    double h = T / M;  // step size
    Eigen::VectorXd yT = y0;  // initial value    
    // Take N RK4 steps:
    for (int step = 0; step < M; step++) {
      // yT is the solution at time t=h*step
      y_next = SystemODE::rk4step(f, h, yT);
      yT = y_next;
    }
    Error(k) = (yT - yT_exact).norm();
    std::cout << std::setw(8) << M << std::setw(20) << Error(k) << std::endl;
  }

  // Estimate convergence rate
  // Get natural logarithm of M by log(M) = log(2)*log2(N).
  Eigen::VectorXd logM =
      std::log(2) * Eigen::VectorXd::LinSpaced(kmax, 1, kmax);
  Eigen::VectorXd coeffs = polyfit(logM, Error.array().log(), 1);
  conv_rate = coeffs(0);

  std::cout << "Convergence rate: "
            << std::round(std::abs(conv_rate)) << std::endl;

  return 0;
}
/* SAM_LISTING_END_0 */
