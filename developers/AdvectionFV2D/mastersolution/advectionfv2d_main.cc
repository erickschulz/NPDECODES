/**
 * @file advectionfv2d_main.cc
 * @brief NPDE homework AdvectionFV2D code
 * @author Philipp Egg
 * @date 01.01.2020
 * @copyright Developed at ETH Zurich
 */

#include <iostream>

#include <Eigen/Core>

#include "advectionfv2d.h"

int main() {
  Eigen::VectorXd v = AdvectionFV2D::dummyFunction(0.0, 0);

  const static Eigen::IOFormat CSVFormat(Eigen::FullPrecision,
                                       Eigen::DontAlignCols, ", ", "\n");
  std::cout << v.transpose().format(CSVFormat) << std::endl;

  return 0;
}
