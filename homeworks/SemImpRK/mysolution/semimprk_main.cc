/**
 * @file semimprk_main.cc
 * @brief NPDE homework SemImpRK code
 * @author Unknown, Oliver Rietmann
 * @date 04.04.2021
 * @copyright Developed at ETH Zurich
 */

#include <Eigen/Core>
#include <iostream>

#include "semimprk.h"

int main() {
  auto f = [](Eigen::Vector3d y) -> Eigen::Vector3d {
    return Eigen::Vector3d(y(0) * y(1), y(1) * y(2), y(2) - y(0));
  };
  auto df = [](Eigen::Vector3d y) {
    Eigen::Matrix3d J;
    J << y(1), y(0), 0.0, 0.0, y(2), y(1), -1.0, 0.0, 1.0;
    return J;
  };
  Eigen::Vector3d y0(1.0, 2.0, 3.0);
  unsigned int M = 10;
  double T = 2.0;

  std::cout << "Test of SolveRosenbrock():"
            << SemImpRK::SolveRosenbrock(f, df, y0, M, T).back().transpose()
            << std::endl;

  double cvgRate = SemImpRK::CvgRosenbrock();
  std::cout << "Convergence rate: " << cvgRate << std::endl;

  return 0;
}
