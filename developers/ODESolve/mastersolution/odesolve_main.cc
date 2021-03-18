/**
 * @file odesolve_main.cc
 * @brief NPDE homework ODESolve code
 * @author ?, Philippe Peter
 * @date 18.03.2021
 * @copyright Developed at ETH Zurich
 */

#include <fstream>
#include <iostream>
#include <vector>

#include "odesolve.h"

int main() {
  auto Psi = [](double h, const double y0) -> double { return (1 + h) * y0; };

  double y0 = 1.0;
  double T = 1.0;

  // Test PsiTilde
  double y1 = ODESolve::PsiTilde(Psi, 1, 0.1, 1.0);
  std::cout << "Test psitilde\n" << y1 << std::endl;

  // Test equidistant integration
  unsigned int M = 8;
  std::cout << "Test equidistant integration" << std::endl;
  std::vector<double> Y1 = ODESolve::OdeIntEqui(Psi, T, y0, M);
  for (int i = 0; i < Y1.size(); ++i) {
    std::cout << Y1[i] << std::endl;
  }
  double rate = ODESolve::TestCvpExtrapolatedEuler();
  std::cout << "\nRate = " << rate << std::endl;

  // Test adaptive integration
  std::cout << "\n\nTest adaptive integration" << std::endl;
  auto [my_t, mysol] =
      ODESolve::OdeIntSsCtrl(Psi, 1, y0, 1.0, 0.01, 10e-5, 10e-5, 10e-5);
  for (int i = 0; i < 8; ++i) {
    std::cout << mysol[i] << std::endl;
  }

  // Test adaptive integration for the Tangent IVP
  auto [tan_t, tan_Y] = ODESolve::SolveTangentIVP();

  std::ofstream solution_file;
  solution_file.open("tangent.csv");

  for (int i = 0; i < tan_Y.size(); ++i) {
    solution_file << tan_t[i] << ", " << tan_Y[i] << "\n";
  }
  solution_file.close();

  std::cout << "Generated " CURRENT_BINARY_DIR "/tangent.csv" << std::endl;
  std::system("python3 " CURRENT_SOURCE_DIR
              "/plot_tangent.py " CURRENT_BINARY_DIR
              "/tangent.csv " CURRENT_BINARY_DIR);

  return 0;
}
