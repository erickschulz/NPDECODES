#include "highresevl.h"

int main() {
  std::cout << "Running driver for discrete evolution in conservation form "
               "with linear reconstruction"
            << std::endl;
  unsigned N = 60; // Number of spatial cells
  double T = 0.5;  // Final time

  // spatial interval for simulation
  double a = -2.5; // left bound for computational interval
  double b = 1.5;  // right bound 

  // Initial state distribution
  auto u0 = [](double x) { return std::sin(x); };
  // Numerical flux function: central flux
  auto F = [](double v, double w) { return 0.25 * (v * v + w * w); };
  // Minmod slope limiter for p.w. linear reconstruction
  // Arguments: mu_{j-1}/h, mu_j/h, mu_{j+1}/h
  auto slopes = [](double a, double b, double c) {
    return std::max(0., std::min(b - a, c - b));
  };
  // Fully discrete evolution 
  ConsFV::highresevl(a, b, N, u0, T, F, slopes);
}
