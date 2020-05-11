/**
 * @file engquistoshernumericalflux_main.cc
 * @brief NPDE homework "EngquistOsherNumericalFlux" code
 * @author Oliver Rietmann
 * @date 23.04.2019
 * @copyright Developed at ETH Zurich
 */

#include <Eigen/Core>
#include <cstdlib>
#include <fstream>
#include <iostream>

#include "engquistoshernumericalflux.h"

/* SAM_LISTING_BEGIN_1 */
const static Eigen::IOFormat CSVFormat(Eigen::FullPrecision,
                                       Eigen::DontAlignCols, ", ", "\n");

double u0(double x) { return (0.0 < x && x <= 1.0) ? 1.0 : -1.0; }

int main() {
  unsigned int N = 100;
  double a = -1.2;
  double b = 2.2;
  double T = 1.0;
  double h = (b - a) / N;

  Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(N, a - 0.5 * h, b - 0.5 * h);
  Eigen::VectorXd uinitial = x.unaryExpr(&u0);
  Eigen::VectorXd ufinal =
      EngquistOsherNumericalFlux::solveCP(a, b, uinitial, T);

#if SOLUTION
  std::ofstream file;
  file.open("ufinal.csv");
  file << x.transpose().format(CSVFormat) << std::endl;
  file << ufinal.transpose().format(CSVFormat) << std::endl;
  file.close();

  std::cout << "Generated " CURRENT_BINARY_DIR "/ufinal.csv" << std::endl;
  std::system("python3 " CURRENT_SOURCE_DIR
              "/plot_solution.py " CURRENT_BINARY_DIR
              "/ufinal.csv " CURRENT_BINARY_DIR "/ufinal.eps");
#else
  //====================
  // Your code goes here
  // Use std::ofstream to write the solution to
  // the file "ufinal.csv". To plot this file
  // you may uncomment the following line:
  // std::system("python3 " CURRENT_SOURCE_DIR "/plot_solution.py "
  // CURRENT_BINARY_DIR "/ufinal.csv " CURRENT_BINARY_DIR "/ufinal.eps");
  //====================
#endif

  return 0;
}
/* SAM_LISTING_END_1 */
