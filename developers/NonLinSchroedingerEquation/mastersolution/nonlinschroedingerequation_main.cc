/**
 * @file nonlinschroedingerequation_main.cc
 * @brief NPDE homework NonLinSchroedingerEquation code
 * @author Oliver Rietmann
 * @date 22.04.2020
 * @copyright Developed at ETH Zurich
 */

#include <iostream>

#include <Eigen/Core>

#include "nonlinschroedingerequation.h"

int main() {
  Eigen::VectorXd v = NonLinSchroedingerEquation::dummyFunction(0.0, 0);

  const static Eigen::IOFormat CSVFormat(Eigen::FullPrecision,
                                       Eigen::DontAlignCols, ", ", "\n");
  std::cout << v.transpose().format(CSVFormat) << std::endl;

  return 0;
}
