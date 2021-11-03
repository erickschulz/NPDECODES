/**
 * @file conslawwithsource_main.cc
 * @brief NPDE exam TEMPLATE CODE FILE
 * @author Oliver Rietmann
 * @date 20.07.2020
 * @copyright Developed at SAM, ETH Zurich
 */

#include <Eigen/Core>
#include <iostream>

#include "conslawwithsource.h"

int main() {
  // Initial data
  auto u0 = [](double x) { return (0.0 <= x && x < 1.0) ? 1.0 : 0.0; };

  // Prepare different spacial resolutions
  const int samples = 8;
  Eigen::VectorXi N(samples);
  N(0) = 10;
  for (int i = 1; i < samples; ++i) N(i) = 2 * N(i - 1);

  // Compute total masses at endtime T=3
  Eigen::VectorXd m3(samples);
  for (int i = 0; i < samples; ++i) {
    Eigen::VectorXd m = ConsLawWithSource::traceMass(u0, N(i));
    m3(i) = m(m.size() - 1);
  }

  // Print N vs. total masses
  Eigen::Matrix<double, samples, 2> table;
  std::cout << "          N        m(3)" << std::endl;
  table.col(0) = N.cast<double>();
  table.col(1) = m3;
  std::cout << table << std::endl;

  return 0;
}
