/**
 * @file systemode_main.cc
 * @brief NPDE homework SystemODE
 * @copyright Developed at ETH Zurich
 */

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
#include <cmath>
#include <iomanip>
#include <iostream>

#include "systemode.h"

/* SAM_LISTING_BEGIN_0 */
int main() {
  // PARAMETERS
  double T = 1;
  int n = 5;

  // INITIAL VALUE
  Eigen::VectorXd y0(2 * n);
  for (int i = 0; i < n; ++i) {
    y0(i) = (i + 1.) / n;
    y0(i + n) = -1;
  }

  // SETUP
  double conv_rate = 0;
  std::cout << std::setw(8) << "M" << std::setw(20) << "Error" << std::endl;

  //====================
  // Your code goes here
  //====================

  std::cout << "Convergence rate: " << std::round(std::abs(conv_rate))
            << std::endl;

  return 0;
}
/* SAM_LISTING_END_0 */
