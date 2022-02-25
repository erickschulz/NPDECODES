/**
 * @file matode_main.cc
 * @brief NPDE homework MatODE code
 * @copyright Developed at ETH Zurich
 */

#include <Eigen/Dense>
#include <iomanip>
#include <iostream>
#include <vector>

#include "matode.h"

int main(int /*argc*/, char** /*argv*/) {
  /* SAM_LISTING_BEGIN_6 */
  double h = 0.01;  // stepsize
  Eigen::Vector3d norms;
  // Build M
  Eigen::Matrix3d M;
  M << 8, 1, 6, 3, 5, 7, 9, 9, 2;
  // Build A
  Eigen::Matrix3d A;
  A << 0, 1, 1, -1, 0, 1, -1, -1, 0;
  Eigen::MatrixXd I = Eigen::Matrix3d::Identity();
  //====================
  // Your code goes here
  //====================
  /* SAM_LISTING_END_6 */
  return 0;
}
