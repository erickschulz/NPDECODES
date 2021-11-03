/**
 * @file gradientflow_main.cc
 * @brief NPDE homework GradientFlow code
 * @author Unknown, Oliver Rietmann
 * @date 04.04.2021
 * @copyright Developed at ETH Zurich
 */

#include <Eigen/Core>
#include <iostream>
#include <vector>

#include "gradientflow.h"

const static Eigen::IOFormat CSVFormat(Eigen::FullPrecision,
                                       Eigen::DontAlignCols, ", ", "\n");

int main() {
  // Parameters and initial condition of the Gradient Flow ODE
  double T = 0.1;
  double lambda = 10.0;
  Eigen::Vector2d d(1.0, 0.0);
  Eigen::Vector2d y0(1.0, 0.0);

  // Approximate exact solution using small timesteps
  int M_ref = 10000;
  std::cout << "T = " << T << ", lambda = " << lambda << std::endl;
  Eigen::Vector2d y_ref =
      GradientFlow::SolveGradientFlow(d, lambda, y0, T, M_ref).back();

  std::cout << "Final value (exact): " << y_ref.transpose().format(CSVFormat)
            << std::endl;

  // Compute error table
  Eigen::VectorXi M(6);
  M << 10, 20, 40, 80, 160, 320;
  std::cout << "Error table:\n";
  std::cout << "M\t error norm" << std::endl;
  for (int i = 0; i < M.size(); ++i) {
    Eigen::Vector2d y_approx =
        GradientFlow::SolveGradientFlow(d, lambda, y0, T, M(i)).back();
    std::cout << M(i) << "\t" << (y_approx - y_ref).norm() << "\t" << std::endl;
  }

  return 0;
}
