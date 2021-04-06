/**
 * @file gradientflow_main.cc
 * @brief NPDE homework GradientFlow code
 * @author Unknown, Oliver Rietmann
 * @date 04.04.2021
 * @copyright Developed at ETH Zurich
 */

#include <Eigen/Core>
#include <cmath>
#include <iostream>
#include <vector>

#include "gradientflow.h"

const static Eigen::IOFormat CSVFormat(Eigen::FullPrecision,
                                       Eigen::DontAlignCols, ", ", "\n");

int main() {
  double T = 0.1;
  double lambda = 10.0;
  int N = 10000;
  std::cout << "T = " << T << ", lambda = " << lambda << std::endl;

  Eigen::Vector2d d(1.0, 0.0);
  Eigen::Vector2d y0(1.0, 0.0);
  std::vector<Eigen::VectorXd> Y =
      GradientFlow::solveGradientFlow(d, lambda, y0, T, N);
  std::cout << "Final value (exact): " << Y.back().transpose().format(CSVFormat)
            << std::endl;

  double exact = Y.back()(0);
  Eigen::VectorXi N_list(6);
  N_list << 10, 20, 40, 80, 160, 320;
  std::cout << "Error table:\n";
  std::cout << "N\terror norm" << std::endl;
  for (int i = 0; i < N_list.size(); ++i) {
    std::vector<Eigen::VectorXd> Y =
        GradientFlow::solveGradientFlow(d, lambda, y0, T, N_list(i));
    double approx = Y.back()(0);
    std::cout << N << "\t" << std::abs(approx - exact) << "\t" << std::endl;
  }

  return 0;
}
