/**
 * @file sspdriver_main.cc
 * @brief NPDE homework ExtendedMUSCL code
 * @author Oliver Rietmann
 * @date 04.08.2019
 * @copyright Developed at ETH Zurich
 */

#include <Eigen/Core>
#include <cmath>
#include <iostream>

#include "extendedmuscl.h"

using namespace ExtendedMUSCL;

int main() {
  // Settings for the ODE
  double T = 1.0;
  double y0 = 1.0;
  double yT_exact = std::exp(T);
  auto f = [](double y) { return y; };

  // Choose time-steps 2^(-4), ..., 2^(-10)
  Eigen::VectorXd tau(7);
  tau << 0x1p-4, 0x1p-5, 0x1p-6, 0x1p-7, 0x1p-8, 0x1p-9, 0x1p-10;

  // Compute error of approx. solution at time T for all timestep-sizes in tau
  int N = tau.size();
  Eigen::VectorXd error(N);
  for (int n = 0; n < N; ++n) {
    int steps = (int)(T / tau(n) + 0.5);
    double y = y0;
    for (int i = 0; i < steps; ++i) y = sspEvolop(f, y, tau(n));
    error(n) = std::abs(yT_exact - y);
  }

  // Print the errors at each timestep
  Eigen::MatrixXd table(3, N);
  table.row(0) = tau;
  table.row(1) = error;
  table.row(2) = error.unaryExpr<double (*)(double)>(&std::log2);
  Eigen::IOFormat tableFormat(2, 0, " ", "\n", " ", " ", " ", " ");
  std::cout << "tau \t error \t log_2(error)" << std::endl;
  std::cout << table.transpose().format(tableFormat) << std::endl;

  return 0;
}
