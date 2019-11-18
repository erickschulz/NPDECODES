/**
 * @file laxwendroffscheme_main.cc
 * @brief NPDE homework "LaxWendroffScheme" code
 * @author Oliver Rietmann
 * @date 29.04.2019
 * @copyright Developed at ETH Zurich
 */

#include "laxwendroffscheme.h"

#include <fstream>

#include <Eigen/Core>

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

  return 0;
}
