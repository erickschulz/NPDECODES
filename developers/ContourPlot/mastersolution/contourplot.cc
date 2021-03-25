#include "contourplot.h"

#include <Eigen/Core>

namespace ContourPlot {

/* SAM_LISTING_BEGIN_0 */
Eigen::Matrix<double, 2, Eigen::Dynamic> crookedEgg() {
#if SOLUTION
  auto gradF = [](Eigen::Vector2d x) -> Eigen::Vector2d {
    return 4.0 * x.squaredNorm() * x - 3.0 * x.cwiseAbs2();
  };
  Eigen::Vector2d y0(1.0, 0.0);
  double T = 4.0;
  return computeIsolinePoints(gradF, y0, T);
#else
  //====================
  // Your code goes here
  // Replace the following dummy return value
  // by the matrix containing the isoline points:
  return Eigen::Matrix<double, 2, 42>::Zero();
  //====================
#endif
}
/* SAM_LISTING_END_0 */

}  // namespace ContourPlot
