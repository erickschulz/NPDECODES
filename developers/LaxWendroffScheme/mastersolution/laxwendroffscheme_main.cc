/**
 * @file laxwendroffscheme_main.cc
 * @brief NPDE homework "LaxWendroffScheme" code
 * @author Oliver Rietmann
 * @date 29.04.2019
 * @copyright Developed at ETH Zurich
 */

#include <fstream>
#include <iostream>

#include <Eigen/Core>

#include "laxwendroffscheme.h"

using namespace LaxWendroffScheme;

const static Eigen::IOFormat CSVFormat(Eigen::FullPrecision,
                                       Eigen::DontAlignCols, ", ", "\n");

int main() {
  Eigen::VectorXi M(6);
  M << 20, 40, 80, 160, 320, 640;

  Eigen::VectorXd error_LaxWendroffRP = numexpLaxWendroffRP(M);
  Eigen::VectorXd error_LaxWendroffSmoothU0 = numexpLaxWendroffSmoothU0(M);
  Eigen::VectorXd error_GodunovSmoothU0 = numexpGodunovSmoothU0(M);

  std::ofstream file;
  file.open("convergence.csv");
  file << M.transpose().format(CSVFormat) << std::endl;
  file << error_LaxWendroffRP.transpose().format(CSVFormat) << std::endl;
  file << error_LaxWendroffSmoothU0.transpose().format(CSVFormat) << std::endl;
  file << error_GodunovSmoothU0.transpose().format(CSVFormat) << std::endl;
  file.close();
  std::cout << "Generated " CURRENT_BINARY_DIR "/convergence.csv" << std::endl;

#if SOLUTION
  std::system("python3 " CURRENT_SOURCE_DIR "/plot.py " CURRENT_BINARY_DIR
              "/convergence.csv " CURRENT_BINARY_DIR "/convergence.eps");
#else
  // To plot from convergence.csv uncomment this:
  // std::system("python3 " CURRENT_SOURCE_DIR "/plot.py " CURRENT_BINARY_DIR
  // "/convergence.csv " CURRENT_BINARY_DIR "/convergence.eps");
#endif

  return 0;
}
