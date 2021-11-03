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

  std::cout << "Convergence MIRK for IVP y' = 1+y^2 " << std::endl;
  // Table header
  std::cout << "M"
            << "\t"
            << "yend"
            << "\t"
            << "err" << std::endl;
  for (unsigned int M = 4; M <= 512; M *= 2) {
    // Solve up to time T = 1, using M equidistant steps
    double yend = MIRK::MIRKSolve(f, df, y0, T, M);
    // Compute error
    double err = std::abs(yex - yend);
    // Print table
    std::cout << M << "\t" << yend << "\t" << err << std::endl;
  }
  return 0;
}
/* SAM_LISTING_END_0 */
