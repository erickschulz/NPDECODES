/**
 * @file exponentialintegrator_main.cc
 * @brief NPDE homework ExponentialIntegrator code
 * @author Unknown, Oliver Rietmann
 * @date 04.04.2021
 * @copyright Developed at ETH Zurich
 */

#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>
/**
 * @file exponentialintegrator_main.cc
 * @brief NPDE homework ExponentialIntegrator code
 * @author Unknown, Oliver Rietmann
 * @date 04.04.2021
 * @copyright Developed at ETH Zurich
 */

#include <Eigen/Core>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

#include "exponentialintegrator.h"

int main() {
  /* SAM_LISTING_BEGIN_0 */
  // Final time
  double T = 1.0;
  // Initial value
  Eigen::VectorXd y0(1);
  y0 << 0.1;
  // Function and Jacobian and exact solution
  auto f = [](const Eigen::VectorXd &y) { return y(0) * (1.0 - y(0)); };
  auto df = [](const Eigen::VectorXd &y) {
    Eigen::VectorXd dfy(1);
    dfy << 1.0 - 2.0 * y(0);
    return dfy;
  };
  double exactyT = y0(0) / (y0(0) + (1.0 - y0(0)) * std::exp(-T));

  // Container for errors
  std::vector<double> error(15);

  // Test many step sizes
  for (int j = 0; j < 15; ++j) {
    int N = std::pow(2, j + 1);
    Eigen::VectorXd y = y0;
    double h = T / N;
    //====================
    // Your code goes here
    // TODO: Perform N timesteps with inital data y0 and store the result in y.
    //====================

    error[j] = std::abs(y(0) - exactyT);
    std::cout << std::left << std::setfill(' ') << std::setw(3)
              << "N = " << std::setw(7) << N << std::setw(8)
              << "Error = " << std::setw(13) << error[j];
    if (j > 0) {
      std::cout << std::left << std::setfill(' ') << std::setw(10)
                << "Approximated order = " << std::log2(error[j - 1] / error[j])
                << std::endl;
    } else
      std::cout << std::endl;
  }
  /* SAM_LISTING_END_0 */
  return 0;
}
