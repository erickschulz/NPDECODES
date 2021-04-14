/**
 * @file sdirk_main.cc
 * @brief NPDE homework SDIRK code
 * @author Unknown, Oliver Rietmann
 * @date 31.03.2021
 * @copyright Developed at ETH Zurich
 */

#include <cstdlib>
#include <iostream>

#include "sdirk.h"

int main() {
  // Compute convergence rates
  double rate = SDIRK::CvgSDIRK();
  std::cout << std::endl << "The rate is " << rate << std::endl;

  // Plot stability domain for gamma = 1.0
  std::system("python3 " CURRENT_SOURCE_DIR
              "/stabdomSDIRK.py " CURRENT_BINARY_DIR "/stabdomSDIRK.eps");
  std::cout << "Generated " CURRENT_BINARY_DIR "/stabdomSDIRK.eps" << std::endl;

  return 0;
}
