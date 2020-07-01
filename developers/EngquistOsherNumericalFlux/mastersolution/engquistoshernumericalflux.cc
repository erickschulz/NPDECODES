/**
 * @file engquistoshernumericalflux.cc
 * @brief NPDE homework "EngquistOsherNumericalFlux" code
 * @author Oliver Rietmann
 * @date 23.04.2019
 * @copyright Developed at ETH Zurich
 */

#include "engquistoshernumericalflux.h"

#include <Eigen/Core>
#include <algorithm>
#include <cmath>

namespace EngquistOsherNumericalFlux {

/* SAM_LISTING_BEGIN_1 */
double EngquistOsherNumFlux(double v, double w) {
  double result;
#if SOLUTION
  auto Fv = [](double v) { return 0 <= v ? std::cosh(v) - 0.5 : 0.5; };
  result = Fv(v) + Fv(-w);
#else
  //====================
  // Your code goes here
  //====================
#endif
  return result;
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
Eigen::VectorXd solveCP(double a, double b, Eigen::VectorXd u0, double T) {
  // Find the maximal speed of propagation
  double A = u0.minCoeff();
  double B = u0.maxCoeff();
  double K = std::max(std::abs(std::sinh(A)), std::abs(std::sinh(B)));
  // Set uniform timestep according to CFL condition
  int N = u0.size();
  double h = (b - a) / N;
  double tau_max = h / K;
  double timesteps = std::ceil(T / tau_max);
  double tau = T / timesteps;

  // Main timestepping loop
#if SOLUTION
  Eigen::VectorXd mu(N);
  for (int i = 0; i < timesteps; ++i) {
    mu.swap(u0);
    double F_minus;
    double F_plus = EngquistOsherNumFlux(mu(0), mu(0));
    for (int j = 0; j < N - 1; ++j) {
      F_minus = F_plus;
      F_plus = EngquistOsherNumFlux(mu(j), mu(j + 1));
      u0(j) = mu(j) - tau / h * (F_plus - F_minus);
    }
    F_minus = F_plus;
    F_plus = EngquistOsherNumFlux(mu(N - 1), mu(N - 1));
    u0(N - 1) = mu(N - 1) - tau / h * (F_plus - F_minus);
  }
#else
  //====================
  // Your code goes here
  //====================
#endif

  return u0;
}
/* SAM_LISTING_END_2 */

}  // namespace EngquistOsherNumericalFlux
