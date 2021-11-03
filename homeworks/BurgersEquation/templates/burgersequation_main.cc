/**
 * @file burgersequation_main.cc
 * @brief NPDE homework BurgersEquation code
 * @author Oliver Rietmann
 * @date 15.04.2019
 * @copyright Developed at ETH Zurich
 */

#include <Eigen/Core>
#include <cstdlib>
#include <fstream>
#include <iostream>

#include "burgersequation.h"

int main() {
  /* SAM_LISTING_BEGIN_1 */
  const unsigned int N = 100;
  Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(N + 1, -1.0, 4.0);
  Eigen::VectorXd mu03 = BurgersEquation::solveBurgersGodunov(0.3, N);
  Eigen::VectorXd mu30 = BurgersEquation::solveBurgersGodunov(3.0, N);

  // Write the solutions to a file that can be used for plotting.
  //====================
  // Your code goes here
  //====================
  /* SAM_LISTING_END_1 */

  /* SAM_LISTING_BEGIN_2 */
  Eigen::Matrix<double, 3, 4> result = BurgersEquation::numexpBurgersGodunov();

  // Write the result to a file that can be used for plotting.
  //====================
  // Your code goes here
  //====================
  /* SAM_LISTING_END_2 */

  return 0;
}
