/**
 * @file rk3prey_main.cc
 * @brief NPDE homework RK3Prey code
 * @copyright Developed at ETH Zurich
 */

#include <Eigen/Dense>
#include <iomanip>
#include <iostream>
#include <limits>
#include "rk3prey.h"

typedef std::numeric_limits<double> dbl;

int main() {
  /* SAM_LISTING_BEGIN_0 */
  // Implementation of butcher scheme
  Eigen::Matrix3d A;
  Eigen::Vector3d b;
  A << 0, 0, 0, 1.0 / 3.0, 0, 0, 0, 2.0 / 3.0, 0;
  b << 0.25, 0, 0.75;

  //====================
  // Your code goes here
  //====================

  // Final time for model
  double T = 10.;

  // Initial value for model
  Eigen::Vector2d y0;
  y0 << 100, 5;

  // Array of number of steps (for convergence study)
  int M[8] = {128, 256, 512, 1024, 2048, 4096, 8192, 16384};

  // Reference "exact" value y(10) at final time T = 10 (approximated)
  Eigen::Vector2d y_ref;
  y_ref << 0.319465882659820, 9.730809352326228;

  // Initialize RK with Butcher scheme
  RK3Prey::RKIntegrator RK(A, b);

  // Start convergence study
  std::cout << std::setw(15) << "M" << std::setw(15) << "error" << std::setw(15)
            << "rate" << std::endl;
  //====================
  // Your code goes here
  //====================
  /* SAM_LISTING_END_0 */
}
