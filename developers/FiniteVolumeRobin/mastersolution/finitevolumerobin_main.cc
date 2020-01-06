/**
 * @file finitevolumerobin_main.cc
 * @brief NPDE homework FiniteVolumeRobin code
 * @author Oliver Rietmann, Erick Schulz
 * @date 01.01.2020
 * @copyright Developed at ETH Zurich
 */

#include <iostream>

#include <Eigen/Core>

#include "finitevolumerobin.h"

int main() {
  Eigen::VectorXd v = FiniteVolumeRobin::dummyFunction(0.0, 0);

  const static Eigen::IOFormat CSVFormat(Eigen::FullPrecision,
                                       Eigen::DontAlignCols, ", ", "\n");
  std::cout << v.transpose().format(CSVFormat) << std::endl;

  return 0;
}
