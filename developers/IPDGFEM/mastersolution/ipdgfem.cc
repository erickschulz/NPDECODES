/**
 * @file ipdgfem.cc
 * @brief NPDE homework IPDGFEM code
 * @author Philippe Peter
 * @date 22.11.2019
 * @copyright Developed at ETH Zurich
 */

#include "ipdgfem.h"

#include <Eigen/Core>

namespace IPDGFEM {

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

}  // namespace IPDGFEM
