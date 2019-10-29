/**
 * @file engquistoshernumericalflux_main.cc
 * @brief NPDE homework "EngquistOsherNumericalFlux" code
 * @author Oliver Rietmann
 * @date 23.04.2019
 * @copyright Developed at ETH Zurich
 */

#include "engquistoshernumericalflux.h"

#include <fstream>
#include <iostream>

#include <Eigen/Core>

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

  /* SOLUTION_BEGIN */
  /* SOLUTION_END */

  return 0;
}
