/**
 * @file newproblem_main.cc
 * @brief NPDE homework NewProblem code
 * @author Oliver Rietmann, Erick Schulz
 * @date 01.01.2020
 * @copyright Developed at ETH Zurich
 */

#include "newproblem.h"

#include <iostream>

#include <Eigen/Core>

const static Eigen::IOFormat CSVFormat(Eigen::FullPrecision,
                                       Eigen::DontAlignCols, ", ", "\n");

int main() {
  Eigen::VectorXd v = NewProblem::dummyFunction(0.0, 0);
  std::cout << v.transpose().format(CSVFormat) << std::endl;

  return 0;
}
