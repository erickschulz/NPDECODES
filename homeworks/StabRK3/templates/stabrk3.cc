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

  //====================
  // Your code goes here
  //====================

  return y;
}
/* SAM_LISTING_END_0 */

/* SAM_LISTING_BEGIN_1 */
void SimulatePredPrey() {
  //====================
  // Your code goes here
  //====================
}

void PrintErrorTable(const Eigen::ArrayXd& M, const Eigen::ArrayXd& error) {
  std::cout << std::setw(15) << "N" << std::setw(15) << "error" << std::setw(15)
            << "rate" << std::endl;
  // Formatted output in C++
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
