/**
 * @file engquistoshernumericalflux.cc
 * @brief NPDE homework "EngquistOsherNumericalFlux" code
 * @author Oliver Rietmann
 * @date 23.04.2019
 * @copyright Developed at ETH Zurich
 */

#include "engquistoshernumericalflux.h"

#include <algorithm>
#include <cmath>

#include <Eigen/Core>

namespace EngquistOsherNumericalFlux {

double EngquistOsherNumFlux(double v, double w) {
  double result;
  /* SOLUTION_BEGIN */
  /* SOLUTION_END */
  return result;
}

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
  /* SOLUTION_BEGIN */
  /* SOLUTION_END */

  return u0;
}

}  // namespace EngquistOsherNumericalFlux
