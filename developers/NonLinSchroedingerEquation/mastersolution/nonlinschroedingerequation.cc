/**
 * @file nonlinschroedingerequation.cc
 * @brief NPDE homework NonLinSchroedingerEquation code
 * @author Oliver Rietmann
 * @date 22.04.2020
 * @copyright Developed at ETH Zurich
 */

#include "nonlinschroedingerequation.h"

#include <Eigen/Core>

namespace NonLinSchroedingerEquation {

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

}  // namespace NonLinSchroedingerEquation
