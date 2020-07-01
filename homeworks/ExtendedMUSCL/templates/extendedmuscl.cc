/**
 * @file extendedmuscl.cc
 * @brief NPDE homework ExtendedMUSCL code
 * @author Oliver Rietmann
 * @date 04.08.2019
 * @copyright Developed at ETH Zurich
 */

#include "extendedmuscl.h"

#include <algorithm>
#include <cmath>
#include <initializer_list>

namespace ExtendedMUSCL {

/* SAM_LISTING_BEGIN_1 */
double logGodunovFlux(double v, double w) {
  double godunov_numerical_flux;
  auto f = [](double u) { return u * (std::log(u) - 1.0); };
  auto df = [](double u) { return std::log(u); };
  //====================
  // Your code goes here
  //====================
  return godunov_numerical_flux;
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_4 */
double limiterMC(double mu_left, double mu_center, double mu_right) {
  double scaled_slope;

  //====================
  // Your code goes here
  //====================

  return scaled_slope;
}
/* SAM_LISTING_END_4 */

}  // namespace ExtendedMUSCL
