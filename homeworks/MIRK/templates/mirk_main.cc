/**
 * @file mirk_main.h
 * @brief NPDE homework MIRK code
 * @author Unknown, Oliver Rietmann
 * @date 04.04.2021
 * @copyright Developed at ETH Zurich
 */

#include <cmath>
#include <iostream>

#include "mirk.h"

/* SAM_LISTING_BEGIN_0 */
int main() {
  // r.h.s
  auto f = [](double y) -> double { return 1 + y * y; };
  // Jacobian of $f$
  auto df = [](double y) -> double { return 2 * y; };
  // Initial data
  const double y0 = 0.;
  // Final time
  const double T = 1.;
  // Exact solution at t = T = 1
  const double yex = tan(T);

  //// PROBLEM h: TEST
  std::cout << "*** PROBLEM h:" << std::endl;
  //====================
  // Your code goes here
  // TODO: problem h: solve IVP y' = f(y) up to T
  //====================
  return 0;
}
/* SAM_LISTING_END_0 */
