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
#if SOLUTION
  if (v >= w)
    godunov_numerical_flux = std::max(f(v), f(w));
  else {
    if (df(v) > 0.0)
      godunov_numerical_flux = f(v);
    else if (df(w) < 0.0)
      godunov_numerical_flux = f(w);
    else
      godunov_numerical_flux = f(1.0);  // note: df(1.0) = 0.0
  }
#else
  //====================
  // Your code goes here
  //====================
#endif
  return godunov_numerical_flux;
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_4 */
double limiterMC(double mu_left, double mu_center, double mu_right) {
  double scaled_slope;

#if SOLUTION
  double sigma_l = 2.0 * (mu_center - mu_left);
  double sigma_c = (mu_right - mu_left) / 2.0;
  double sigma_r = 2.0 * (mu_right - mu_center);

  std::initializer_list<double> sigma_list{sigma_l, sigma_c, sigma_r};

  double min = std::min(sigma_list);
  if (min > 0.0)
    scaled_slope = min;
  else {
    double max = std::max(sigma_list);
    if (max < 0.0)
      scaled_slope = max;
    else
      scaled_slope = 0.0;
  }
#else
  //====================
  // Your code goes here
  //====================
#endif

  return scaled_slope;
}
/* SAM_LISTING_END_4 */

}  // namespace ExtendedMUSCL
