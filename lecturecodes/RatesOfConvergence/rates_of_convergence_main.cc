/** @file
 * @brief functions for estimating order of convergence
 * @author Ralf Hiptmair
 * @date March 2019
 */

#include <iostream>
#include "estimated_order_convergence.h"

int main(int /*argc*/, char** /*argv*/) {
  Eigen::VectorXd N(7);
  N << 2, 4, 8, 16, 32, 64, 128;

  Eigen::VectorXd err(7);
  err << 0.5, 0.26, 0.123, 0.060, 0.032, 0.014, 0.00711;

  double rate = rate_of_convergence::eoc(N, err);
  std::cout << "rate extracted from data = " << rate << std::endl;

  return 0;
}
