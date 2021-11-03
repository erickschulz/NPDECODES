/**
 * @file numpdesetup.cc
 * @brief NPDE homework NumPDESetup code
 * @author Oliver Rietmann, Erick Schulz
 * @date 01.01.2020
 * @copyright Developed at ETH Zurich
 */

#include "numpdesetup.h"

#include <Eigen/Core>
#include <iostream>

namespace NumPDESetup {

/* SAM_LISTING_BEGIN_1 */
Eigen::VectorXd dummyFunction(double x, int n) {
  // Appears only in mastersolution
  std::cout << "NumPDESetup: master solution code" << std::endl;
  return Eigen::VectorXd::Constant(n, x);
}
/* SAM_LISTING_END_1 */

}  // namespace NumPDESetup
