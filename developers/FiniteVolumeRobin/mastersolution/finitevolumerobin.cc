/**
 * @file finitevolumerobin.cc
 * @brief NPDE homework FiniteVolumeRobin code
 * @author Oliver Rietmann, Erick Schulz
 * @date 01.01.2020
 * @copyright Developed at ETH Zurich
 */

#include "finitevolumerobin.h"

#include <Eigen/Core>

namespace FiniteVolumeRobin {

/* SAM_LISTING_BEGIN_1 */
Eigen::VectorXd dummyFunction(double x, int n) {
#if SOLUTION
  // Appears only in mastersolution
  return Eigen::Vector2d::Constant(1.0);
#else
  // Appears only in mysolution and templates
  return Eigen::Vector2d::Zero();
#endif
}
/* SAM_LISTING_END_1 */

}  // namespace FiniteVolumeRobin
