/**
 * @file initcondlv_main.cc
 * @brief NPDE homework InitCondLV code
 * @author lfilippo, tille, jgacon, dcasati
 * @copyright Developed at ETH Zurich
 */

#include <Eigen/Dense>
#include <iostream>
#include <utility>

#include "initcondlv.h"

int main() {
  /* Compute initial conditions u0, v0 such that solution has period T=5
   * using inital approximation [3,2]^T */
  /* SAM_LISTING_BEGIN_2 */
  Eigen::Vector2d y(3, 2);  // Initial guess
  double T = 5;             // Period

  //====================
  // Your code goes here
  //====================

  // Calculate the evolution of the Lotka-Volterra ODE up to
  // time 100, using the initial value y, found by the Newton iterations above.
  // If the evolution is indeed 5-periodic, then the solution at time 100
  // should have the value y. Print a warning message if this is not the case.
  //====================
  // Your code goes here
  //====================

  // Testing obtained solution against reference "exact" solution
  Eigen::Vector2d ref(3.1098751029156, 2.08097564048345);
  double tol = 1.0e-8;
  double error = (ref - y).lpNorm<Eigen::Infinity>();
  if (std::abs(error) > tol) {
    std::cout
        << " Error w.r.t. to reference solution is GREATER than tol=1.0e-8."
        << std::endl;
  } else {
    std::cout
        << " Error w.r.t. to reference solution is SMALLER than tol=1.0e-8."
        << std::endl;
  }

  /* SAM_LISTING_END_2 */
  return 0;
}
