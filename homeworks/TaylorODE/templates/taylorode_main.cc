#include <iostream>

#include "taylorode.h"

/**
 * @file taylorode_main.cc
 * @brief NPDE homework TaylorODE
 * @author Philippe Peter
 * @date 24.03.2021
 * @copyright Developed at ETH Zurich
 */

int main() {
  // Verify rate of convergence
  double convRate = TaylorODE::TestCvgTaylorMethod();
  std::cout << "Estimated rate of convergence: " << convRate << std::endl;
  std::cout << "Expected rate of convergence: " << 3.0 << std::endl;
  return 0;
}
