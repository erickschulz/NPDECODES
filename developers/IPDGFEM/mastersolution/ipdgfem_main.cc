/**
 * @file ipdgfem_main.cc
 * @brief NPDE homework IPDGFEM code
 * @author Philippe Peter
 * @date 22.11.2019
 * @copyright Developed at ETH Zurich
 */

#include "ipdgfem.h"

#include <iostream>

#include <Eigen/Core>

const static Eigen::IOFormat CSVFormat(Eigen::FullPrecision,
                                       Eigen::DontAlignCols, ", ", "\n");

int main() {
  Eigen::VectorXd v = IPDGFEM::dummyFunction(0.0, 0);
  std::cout << v.transpose().format(CSVFormat) << std::endl;

  return 0;
}
