/**
 * @file stabrk3.cc
 * @brief NPDE homework StabRK3 code
 * @author Unknown, Oliver Rietmann, Philippe Peter
 * @date 13.04.2021
 * @copyright Developed at ETH Zurich
 */

#include "stabrk3.h"

#include <Eigen/Core>

#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

namespace StabRK3 {

/* SAM_LISTING_BEGIN_0 */
Eigen::Vector2d PredPrey(Eigen::Vector2d y0, double T, unsigned int M) {
  double h = T / M;
  Eigen::Vector2d y = y0;

  auto f = [](Eigen::Vector2d y) -> Eigen::Vector2d {
    return {(1 - y(1)) * y(0), (y(0) - 1) * y(1)};
  };

  for (int j = 0; j < M; ++j) {
    Eigen::Vector2d k1 = f(y);
    Eigen::Vector2d k2 = f(y + h * k1);
    Eigen::Vector2d k3 = f(y + (h / 4.) * k1 + (h / 4.) * k2);
    y = y + (h / 6.) * k1 + (h / 6.) * k2 + (2. * h / 3.) * k3;
  }

  return y;
}
/* SAM_LISTING_END_0 */

/* SAM_LISTING_BEGIN_1 */
void SimulatePredPrey() {
  // Parameters
  double T = 1.0;
  Eigen::Vector2d y0(100.0, 1.0);

  //(Approximate) reference solution
  Eigen::Vector2d y_ref = PredPrey(y0, T, std::pow(2, 14));

  Eigen::ArrayXd error(12);
  Eigen::ArrayXd M(12);
  M << 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192;

  // Compute errors
  for (int i = 0; i < M.size(); ++i) {
    Eigen::Vector2d y = PredPrey(y0, T, M(i));
    error(i) = (y - y_ref).norm();
  }

  // Print error table
  PrintErrorTable(M, error);
}

void PrintErrorTable(const Eigen::ArrayXd& M, const Eigen::ArrayXd& error) {
  std::cout << std::setw(15) << "N" << std::setw(15) << "error" << std::setw(15)
            << "rate" << std::endl;

  for (unsigned int i = 0; i < M.size(); ++i) {
    std::cout << std::setw(15) << M(i) << std::setw(15) << error(i);
    if (i > 0) {
      std::cout << std::setw(15) << std::log2(error(i - 1) / error(i));
    }
    std::cout << std::endl;
  }
}
/* SAM_LISTING_END_1 */

}  // namespace StabRK3
