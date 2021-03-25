/**
 * @file contourplot.cc
 * @brief NPDE homework ContourPlot code
 * @author Unknown, Oliver Rietmann
 * @date 25.03.2021
 * @copyright Developed at ETH Zurich
 */

#include "contourplot.h"

#include <Eigen/Core>

namespace ContourPlot {

/* SAM_LISTING_BEGIN_0 */
Eigen::Matrix<double, 2, Eigen::Dynamic> crookedEgg() {
  //====================
  // Your code goes here
  // Replace the following dummy return value
  // by the matrix containing the isoline points:
  return Eigen::Matrix<double, 2, 42>::Zero();
  //====================
}
/* SAM_LISTING_END_0 */

}  // namespace ContourPlot
