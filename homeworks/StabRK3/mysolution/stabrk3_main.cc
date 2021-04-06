/**
 * @file stabrk3_main.cc
 * @brief NPDE homework StabRK3 code
 * @author Oliver Rietmann
 * @date 04.04.2021
 * @copyright Developed at ETH Zurich
 */

#include <Eigen/Core>
#include <iomanip>
#include <iostream>
#include <vector>

#include "stabrk3.h"

const static Eigen::IOFormat CSVFormat(Eigen::FullPrecision,
                                       Eigen::DontAlignCols, ", ", "\n");

int main() {
  double T = 1.0;
  Eigen::Vector2d y0(100.0, 1.0);
  Eigen::Vector2d yT_reference = StabRK3::predPrey(y0, T, 16384);
  std::cout << "Solution Computed by predPrey(): "
            << yT_reference.transpose().format(CSVFormat) << std::endl;

  // Vector of number of steps (for convergence study)
  std::vector<unsigned int> N_list = {4,   8,   16,   32,   64,   128,
                                      256, 512, 1024, 2048, 4096, 8192};

  // Compute approximations and take last one as reference
  std::vector<Eigen::Vector2d> yT_list = StabRK3::simulatePredPrey(N_list);

  // Compute the error table
  std::cout << std::setw(15) << "N" << std::setw(15) << "error" << std::setw(15)
            << "rate" << std::endl;
  double error_old;
  for (unsigned int j = 0; j < N_list.size(); ++j) {
    double error = (yT_list[j] - yT_reference).norm();
    std::cout << std::setw(15) << N_list[j] << std::setw(15) << error;
    if (j > 0) {
      std::cout << std::setw(15) << log2(error_old / error);
    }
    error_old = error;
    std::cout << std::endl;
  }

  return 0;
}
