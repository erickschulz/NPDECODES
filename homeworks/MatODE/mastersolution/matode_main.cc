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
  // Build Q
  Eigen::HouseholderQR<Eigen::Matrix3d> qr(3, 3);
  qr.compute(M);
  Eigen::MatrixXd Q = qr.householderQ();
  // Set initial conditions
  Eigen::Matrix3d Meeul = Q, Mieul = Q, Mimp = Q;

  std::vector<int> sep = {5, 15};
  std::cout << std::setw(sep[0]) << "step" << std::setw(sep[1]) << "exp. Eul"
            << std::setw(sep[1]) << "imp. Eul" << std::setw(sep[1]) << "Mid-Pt"
            << std::endl;
  // Norm of Y'Y-I for initial value
  std::cout << std::setw(sep[0]) << "0" << std::setw(sep[1])
            << (Meeul.transpose() * Meeul - I).norm() << std::setw(sep[1])
            << (Mieul.transpose() * Mieul - I).norm() << std::setw(sep[1])
            << (Mimp.transpose() * Mimp - I).norm() << std::endl;
  // Norm of Y'Y-I for 20 steps
  for (unsigned int j = 0; j < 20; ++j) {
    Meeul = MatODE::eeulstep(A, Meeul, h);
    Mieul = MatODE::ieulstep(A, Mieul, h);
    Mimp = MatODE::impstep(A, Mimp, h);

    norms[0] = (Meeul.transpose() * Meeul - I).norm();
    norms[1] = (Mieul.transpose() * Mieul - I).norm();
    norms[2] = (Mimp.transpose() * Mimp - I).norm();

    std::cout << std::setw(sep[0]) << j + 1 << std::setw(sep[1]) << norms[0]
              << std::setw(sep[1]) << norms[1] << std::setw(sep[1]) << norms[2]
              << std::endl;
  }
  /* SAM_LISTING_END_6 */
  return 0;
}
