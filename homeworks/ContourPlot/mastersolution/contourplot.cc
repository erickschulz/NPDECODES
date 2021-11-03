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
  auto gradF = [](Eigen::Vector2d x) -> Eigen::Vector2d {
    return 4.0 * x.squaredNorm() * x - 3.0 * x.cwiseAbs2();
  };
  Eigen::Vector2d y0(1.0, 0.0);
  double T = 4.0;
  return computeIsolinePoints(gradF, y0, T);
}
/* SAM_LISTING_END_0 */

}  // namespace ContourPlot
