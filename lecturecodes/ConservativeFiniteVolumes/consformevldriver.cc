// * Demonstration code for course Numerical Methods for Partial Differential
// Equations
// * Author: R. Hiptmair, SAM, ETH Zurich
// * Date: May 2020

#include "consformevl.h"
#include "numexp_runner.h"

int main(int /*argc*/, char ** /*argv*/) {
  std::cout << "Model problem: 1D Burgers equation" << std::endl;
  // Spatial computational domain. Note that the initial data are confined to
  // [0,1]. Therefore the maximal speed of propagation "to the right" will be 1.
  // This means that the exact solution  will always be supported in [0,5] for
  // times in [0,4].
  const double _a = -1.0;
  const double _b = 5.0;
  // Final time for simulation
  const double _T = 4.0;

  {
    auto evl = [&](double a, double b, double N, double T) -> Eigen::VectorXd {
      return ConsFV::consformevl(a, b, N, box, T, nfn_lf_burger);
    };
    consform_compute(evl, "burgers_box_rusanov.csv", _T, _a, _b);
  }
  {
    auto evl = [&](double a, double b, double N, double T) -> Eigen::VectorXd {
      return ConsFV::consformevl(a, b, N, box, T, nfn_god_burger);
    };
    consform_compute(evl, "burgers_box_godunov.csv", _T, _a, _b);
  }
  {
    auto evl = [&](double a, double b, double N, double T) -> Eigen::VectorXd {
      return ConsFV::consformevl(a, b, N, wedge, T, nfn_lf_burger);
    };
    consform_compute(evl, "burgers_wedge_rusanov.csv", _T, _a, _b);
  }
  {
    auto evl = [&](double a, double b, double N, double T) -> Eigen::VectorXd {
      return ConsFV::consformevl(a, b, N, wedge, T, nfn_god_burger);
    };
    consform_compute(evl, "burgers_wedge_godunov.csv", _T, _a, _b);
  }
  {
    auto evl = [&](double a, double b, double N, double T) -> Eigen::VectorXd {
      return ConsFV::consformevl(a, b, N, bump, T, nfn_lf_burger);
    };
    consform_compute(evl, "burgers_bump_rusanov.csv", _T, _a, _b);
  }
  {
    auto evl = [&](double a, double b, double N, double T) -> Eigen::VectorXd {
      return ConsFV::consformevl(a, b, N, bump, T, nfn_god_burger);
    };
    consform_compute(evl, "burgers_bump_godunov.csv", _T, _a, _b);
  }
}
