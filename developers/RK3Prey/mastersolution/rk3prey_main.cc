#include <Eigen/Dense>
#include <iomanip>
#include <iostream>
#include <vector>

#include "rkintegrator.h"

//! \file rk3prey.cc Solve prey/predator model with RK-SSM method

int main() {
  /* SAM_LISTING_BEGIN_0 */
  // Implementation of butcher scheme
  unsigned int s = 3;
  Eigen::MatrixXd A(s, s);
  Eigen::VectorXd b(s);
  A << 0, 0, 0, 1. / 3., 0, 0, 0, 2. / 3., 0;
  b << 1. / 4., 0, 3. / 4.;

  // Coefficients and handle for prey/predator model
  double alpha1 = 3.;
  double alpha2 = 2.;
  double beta1 = 0.1;
  double beta2 = 0.1;
#if SOLUTION
  auto f = [&alpha1, &alpha2, &beta1, &beta2](const Eigen::VectorXd& y) {
    auto temp = y;
    temp(0) *= alpha1 - beta1 * y(1);
    temp(1) *= -alpha2 + beta2 * y(0);
    return temp;
  };
#else   // TEMPLATE
  // TODO: implement functor $f$ for rhs of $y(t)' = f(y(t))$
#endif  // TEMPLATE

  // Dimension of state space
  unsigned int d = 2;

  // Final time for model
  double T = 10.;

  // Initial value for model
  Eigen::VectorXd y0(d);
  y0 << 100, 5;

  // Array of number of steps (for convergence study)
  std::vector<unsigned int> N = {128, 256, 512, 1024, 2048, 4096, 8192, 16384};

  // Exact value y(10) at final time T = 10 (approximated)
  Eigen::VectorXd yex(d);
  yex << 0.319465882659820, 9.730809352326228;

  // Initialize RK with Butcher scheme
  RK3Prey::RKIntegrator<Eigen::VectorXd> RK(A, b);

  // Start convergence study
  std::cout << std::setw(15) << "N" << std::setw(15) << "error" << std::setw(15)
            << "rate" << std::endl;
#if SOLUTION
  double errold;
  for (unsigned int i = 0; i < N.size(); ++i) {
    auto res = RK.solve(f, T, y0, N[i]);
    double err = (res.back() - yex).norm();
    std::cout << std::setw(15) << N[i] << std::setw(15) << err;
    if (i > 0) {
      std::cout << std::setw(15) << log2(errold / err);
    }
    errold = err;
    std::cout << std::endl;
  }
#else   // TEMPLATE
  // TODO: tabulate error for each $N$
#endif  // TEMPLATE
  /* SAM_LISTING_END_0 */
}
