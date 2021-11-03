/**
 * @file numpdesetup_main.cc
 * @brief NPDE homework NumPDESetup code
 * @author Ralf Hiptmair
 * @date 17.02.2021
 * @copyright Developed at ETH Zurich
 */

#include <Eigen/Core>
#include <iostream>

#include "numpdesetup.h"

/* SAM_LISTING_BEGIN_2 */
int main(int /*argc*/, char** /*argv*/) {
  std::cout << "Dummy \"homework problem\" for NumPDE course" << std::endl;
  Eigen::VectorXd v = NumPDESetup::dummyFunction(42.0, 10);

  const static Eigen::IOFormat CSVFormat(Eigen::FullPrecision,
                                         Eigen::DontAlignCols, ", ", "\n");
  std::cout << v.transpose().format(CSVFormat) << std::endl;
  return 0;
}
/* SAM_LISTING_END_2 */
