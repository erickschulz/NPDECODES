#include "consformevl.h"

#include <cmath>
#include <iostream>

int main(int /*argc*/, char **/*argv*/) {
  std::cout << "Running driver for discrete evolution in conservation form"
            << std::endl;
  // Test script for discrete evolution in conservation form
  double a = -2 * M_PI;
  double b = 2 * M_PI;
  unsigned N = 100;
  auto u0 = [](double x) { return std::sin(x); };
  double T = 100.;
  auto F = [](double v, double w) { return 0.25 * (v * v + w * w); };

  ConsFV::consformevl(a, b, N, u0, T, F);
}
