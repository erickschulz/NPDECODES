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
#if SOLUTION
  // Table header
  std::cout << "N"
            << "\t"
            << "yend"
            << "\t"
            << "err" << std::endl;
  for (unsigned int N = 4; N < 512; N = N << 1) {
    // Solve up to time T = 1, using N equidistant steps
    double yend = MIRK::mirkSolve(f, df, y0, T, N);
    // Compute error
    double err = std::abs(yex - yend);
    // Print table
    std::cout << N << "\t" << yend << "\t" << err << std::endl;
  }
#else
  //====================
  // Your code goes here
  // TODO: problem h: solve IVP y' = f(y) up to T
  //====================
#endif
  return 0;
}
/* SAM_LISTING_END_0 */
