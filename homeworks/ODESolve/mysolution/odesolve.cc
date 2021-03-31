/**
 * @file odesolve.cc
 * @brief NPDE homework ODESolve code
 * @author ?, Philippe Peter
 * @date 18.03.2021
 * @copyright Developed at ETH Zurich
 */

#include "odesolve.h"

#include <Eigen/Core>
#include <cmath>
#include <iostream>
#include <vector>

#include "polyfit.h"

namespace ODESolve {

/* SAM_LISTING_BEGIN_2 */
double TestCvpExtrapolatedEuler() {
  double conv_rate;
  double T = 1.0;
  double y0 = 0.0;
  auto f = [](double y) -> double { return 1.0 + y * y; };

  // TODO : tabulate the values of the error corresponding to
  // \tilde{\psi}, where \psi is the explicit Euler method.
  // return the empirical convergence rate using polyfit.
  // Hint: first define a lambda for \psi. Then use psitilde to obtain a
  // suitable input for odeintequi.

  // ===================
  // Your code goes here
  // ===================
  return conv_rate;
}

/* SAM_LISTING_BEGIN_4 */
std::pair<std::vector<double>, std::vector<double>> SolveTangentIVP() {
  auto f = [](double y) -> double { return 1.0 + y * y; };
  double y0 = 0.0;

  // TODO: run the adaptive integration algorithm
  // ===================
  // Your code goes here
  // ===================
  // dummy vectors (if size 0, the plotting script crashes)
  std::vector<double> t(10, 0.0);
  std::vector<double> Y(10, 1.0);
  return {t, Y};
}
/* SAM_LISTING_END_4 */

}  // namespace ODESolve