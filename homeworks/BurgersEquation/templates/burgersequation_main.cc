/**
 * @file burgersequation_main.cc
 * @brief NPDE homework BurgersEquation code
 * @author Oliver Rietmann
 * @date 15.04.2019
 * @copyright Developed at ETH Zurich
 */

#include "burgersequation.h"

#include <fstream>
#include <iostream>

#include <Eigen/Core>

int main() {
  const unsigned int N = 100;
  Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(N + 1, -1.0, 4.0);
  Eigen::VectorXd mu03 = BurgersEquation::solveBurgersGodunov(0.3, N);
  Eigen::VectorXd mu30 = BurgersEquation::solveBurgersGodunov(3.0, N);

  // Write the solutions to a file that can be used for plotting.
  /* SOLUTION_BEGIN */
  /* SOLUTION_END */

  Eigen::Matrix<double, 3, 4> result = BurgersEquation::numexpBurgersGodunov();

  // Write the result to a file that can be used for plotting.
  /* SOLUTION_BEGIN */
  /* SOLUTION_END */

  return 0;
}
