/**
 * @file numpdesetup.cc
 * @brief NPDE homework NumPDESetup code
 * @author Oliver Rietmann, Erick Schulz
 * @date 01.01.2020
 * @copyright Developed at ETH Zurich
 */

#include "numpdesetup.h"

#include <Eigen/Core>

namespace NumPDESetup {

/* SAM_LISTING_BEGIN_1 */
Eigen::VectorXd dummyFunction(double x, int n) {
#if SOLUTION
  // Appears only in mastersolution
  return Eigen::VectorXd::Constant(n,x);
#else
  // Appears only in mysolution and templates
  return Eigen::Vector2d::Zero();
#endif
}
/* SAM_LISTING_END_1 */

}  // namespace NumPDESetup
