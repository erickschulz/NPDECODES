/**
 * @file semimprk.cc
 * @brief NPDE homework SemImpRK code
 * @author Unknown, Oliver Rietmann
 * @date 04.04.2021
 * @copyright Developed at ETH Zurich
 */

#include "semimprk.h"

#include <Eigen/Core>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

#include "../../../lecturecodes/helperfiles/polyfit.h"

namespace SemImpRK {

/* SAM_LISTING_BEGIN_0 */
double CvgRosenbrock() {
  double cvgRate = 0.0;
  // Use polyfit to estimate the rate of convergence
  // for SolveRosenbrock.
  //====================
  // Your code goes here
  //====================
  return cvgRate;
}
/* SAM_LISTING_END_0 */

}  // namespace SemImpRK
