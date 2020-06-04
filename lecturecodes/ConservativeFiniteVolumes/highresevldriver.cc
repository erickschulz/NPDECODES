// Demonstration code for course Numerical Methods for Partial Differential
// Equations Author: R. Hiptmair, SAM, ETH Zurich Date: May 2020

#include "highresevl.h"
#include "numexp_runner.h"

int main() {
  std::cout << "Model problem: 1D Burgers equation" << std::endl;
  std::cout
      << "Numerical simulation based on high-resolution FV with slope limiting"
      << std::endl;

  // Spatial computational domain. Note that the initial data are confined to
  // [0,1]. Therefore the maximal speed of propagation "to the right" will be 1.
  // This means that the exact solution  will always be supported in [0,5] for
  // times in [0,4].
  const double _a = -1.0;
  const double _b = 5.0;
  // Final time for simulation
  const double _T = 4.0;

  // No linear reconstruction: for testing purposes
  auto zeroslope = [](double a, double b, double c) { return 0.0; };
  
  // Minmod slope limiter for p.w. linear reconstruction
  // Arguments: mu_{j-1}/h, mu_j/h, mu_{j+1}/h
  auto minmod = [](double a, double b, double c) {
    return std::max(0., std::min(b - a, c - b));
  };
  
  {
    auto evl = [&](double a, double b, double N, double T) -> Eigen::VectorXd {
      return ConsFV::highresevl(a, b, N, box, T, nfn_lf_burger, minmod);
    };
    consform_compute(evl, "burgers_hr_rusanov.csv", _T, _a, _b);
  }
}
