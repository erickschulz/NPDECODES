/**
 * @file odesolve.cc
 * @brief NPDE homework ODESolve code
 * @author ?, Philippe Peter
 * @date 18.03.2021
 * @copyright Developed at ETH Zurich
 */

#include "odesolve.h"

#include <Eigen/Core>
#include <cmath>
#include <iostream>
#include <vector>

#include "../../../lecturecodes/helperfiles/polyfit.h"

namespace ODESolve {

/* SAM_LISTING_BEGIN_2 */
double TestCvpExtrapolatedEuler() {
  double conv_rate;
  double T = 1.0;
  double y0 = 0.0;
  auto f = [](double y) -> double { return 1.0 + y * y; };

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
  Eigen::ArrayXd M(11);

  std::cout << "Error table for equidistant steps:" << std::endl;
  std::cout << "M"
            << "\t"
            << "Error" << std::endl;

  for (int i = 0; i < 11; ++i) {
    M(i) = std::pow(2, i + 2);
    double yT = OdeIntEqui(Psi_tilde, T, y0, M(i)).back();
    err(i) = std::abs(yT - y_ex1);
    std::cout << M(i) << "\t" << err(i) << std::endl;
  }
  // compute fitted rate
  Eigen::VectorXd coeffs = polyfit(M.log(), err.log(), 1);
  conv_rate = -coeffs(0);
  return conv_rate;
}

/* SAM_LISTING_BEGIN_4 */
std::pair<std::vector<double>, std::vector<double>> SolveTangentIVP() {
  auto f = [](double y) -> double { return 1.0 + y * y; };
  double y0 = 0.0;

  double T = 1.5;
  unsigned p = 1;
  double h0 = 1. / 100.;

  auto Psi = [&f](double h, double y0) -> double { return y0 + h * f(y0); };
  // run adaptive algoritm
  return OdeIntSsCtrl(Psi, p, y0, T, h0, 10e-4, 10e-6, 10e-5);

}
/* SAM_LISTING_END_4 */

}  // namespace ODESolve
