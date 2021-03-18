/**
 * @file odesolve.cc
 * @brief NPDE homework ODESolve code
 * @author ?, Philippe Peter
 * @date 18.03.2021
 * @copyright Developed at ETH Zurich
 */

#include <Eigen/Core>

#include <cmath>
#include <iostream>
#include <vector>

#include "odesolve.h"
#include "polyfit.h"

namespace ODESolve {

/* SAM_LISTING_BEGIN_2 */
double TestCvpExtrapolatedEuler() {
  double conv_rate;
  double T = 1.0;
  double y0 = 0.0;
  auto f = [](double y) -> double { return 1.0 + y * y; };

#if SOLUTION
  auto Psi = [&f](double h, double y0) -> double { return y0 + h * f(y0); };
  unsigned p = 1;

  // lambda corresponding to \tilde{\psi}
  auto Psi_tilde = [&Psi, &f, &p](double h, double y0) -> double {
    return PsiTilde(Psi, p, h, y0);
  };

  // exact value
  double y_ex1 = std::tan(T);

  // values for convergence study
  Eigen::ArrayXd err(11);
  Eigen::ArrayXd N(11);

  std::cout << "Error table for equidistant steps:" << std::endl;
  std::cout << "N"
            << "\t"
            << "Error" << std::endl;

  for (int i = 0; i < 11; ++i) {
    N(i) = std::pow(2, i + 2);
    double yT = OdeIntEqui(Psi_tilde, T, y0, N(i)).back();
    err(i) = std::abs(yT - y_ex1);
    std::cout << N(i) << "\t" << err(i) << std::endl;
  }
  // compute fitted rate
  Eigen::VectorXd coeffs = polyfit(N.log(), err.log(), 1);
  conv_rate = -coeffs(0);
#else
  // TODO : tabulate the values of the error corresponding to
  // \tilde{\psi}, where \psi is the explicit Euler method.
  // return the empirical convergence rate using polyfit.
  // Hint: first define a lambda for \psi. Then use psitilde to obtain a
  // suitable input for odeintequi.

  // ===================
  // Your code goes here
  // ===================
#endif
  return conv_rate;
}

/* SAM_LISTING_BEGIN_4 */
std::pair<std::vector<double>, std::vector<double>> SolveTangentIVP() {
  auto f = [](double y) -> double { return 1.0 + y * y; };
  double y0 = 0.0;

#if SOLUTION
  double T = 1.5;
  unsigned p = 1;
  double h0 = 1. / 100.;

  auto Psi = [&f](double h, double y0) -> double { return y0 + h * f(y0); };
  // run adaptive algoritm
  return OdeIntSsCtrl(Psi, p, y0, T, h0, 10e-4, 10e-6, 10e-5);

#else
  // TODO: run the adaptive integration algorithm
  // ===================
  // Your code goes here
  // ===================
  // dummy vectors (if size 0, the plotting script crashes)
  std::vector<double> t(10, 0.0);
  std::vector<double> Y(10, 1.0);
  return {t, Y};
#endif
}
/* SAM_LISTING_END_4 */

}  // namespace ODESolve