/**
 * @file matode_main.cc
 * @brief NPDE homework MatODE code
 * @copyright Developed at ETH Zurich
 */

#include <Eigen/Dense>
#include <iomanip>
#include <iostream>

#include "matode.h"

int main() {
  /* SAM_LISTING_BEGIN_6 */
  double h = 0.01;  // stepsize
  Eigen::Vector3d norms;

  // Build M
  Eigen::Matrix3d M;
  M << 8, 1, 6, 3, 5, 7, 9, 9, 2;

  //====================
  // Your code goes here
  //====================

  /* SAM_LISTING_END_6 */
}
