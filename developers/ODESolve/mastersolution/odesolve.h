#ifndef ODESOLVE_H
#define ODESOLVE_H

/**
 * @file odesolve.h
 * @brief NPDE homework ODESolve code
 * @author ?, Philippe Peter
 * @date 18.03.2021
 * @copyright Developed at ETH Zurich
 */

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

namespace ODESolve {

/**
 * @brief Evolves the vector y0 using the evolution operator \tilde{\Psi} with
 * step-size y0.
 * @tparam DiscEvlOp type for evolution operator (e.g. lambda
 * function type)
 * @param Psi original evolution operator
 * @param p parameter p for construction of Psi tilde
 * @param h step-size
 * @param y0 previous step
 * @return Evolved step \f$ \tilde{\Psi}^h y0 \f$
 */
/* SAM_LISTING_BEGIN_0 */
template <class DiscEvlOp>
double PsiTilde(const DiscEvlOp& Psi, unsigned int p, double h, double y0) {
#if SOLUTION
  return (Psi(h, y0) - std::pow(2, p) * Psi(h / 2., Psi(h / 2., y0))) /
         (1. - std::pow(2, p));
#else
  // TODO: implement psi tilde operator and replace the dummy return value y0
  //====================
  // Your code goes here
  //====================
  return y0;
#endif
}
/* SAM_LISTING_END_0 */

/**
 * @brief Evolves the state y0 using the evolution operator from time 0 to T
 * using equidistant steps.
 * @tparam DiscEvlOp type for evolution operator (e.g. lambda function type)
 * @param Psi original evolution operator
 * @param T final time
 * @param y0 initial data
 * @param M number of steps
 * @return vector of data at equidistant time points.
 */
/* SAM_LISTING_BEGIN_1 */
template <class DiscEvlOp>
std::vector<double> OdeIntEqui(const DiscEvlOp& Psi, double T, double y0,
                               int M) {
  double h = T / M;
  double t = 0.;
  std::vector<double> Y;

#if SOLUTION
  Y.reserve(M + 1);
  Y.push_back(y0);
  double y = y0;

  while (t < T) {
    y = Psi(h, y);
    Y.push_back(y);
    t += std::min(T - t, h);
  }
#else
  // TODO: implement equidistant solver for evolution operaotr Psi
  //====================
  // Your code goes here
  //====================
  Y.resize(M);
#endif

  return Y;
}
/* SAM_LISTING_END_1 */

double TestCvpExtrapolatedEuler();

/**
 * @brief Evolves the vector y0 using the evolution operator from time 0 to T
 * using adaptive error control
 * @tparam DiscEvlOp type of the evolution operator
 * @param Psi original evolution operator
 * @param p parameter p for construction of Psi tilde
 * @param y0 initial data
 * @param T final time
 * @param h0 initial step size
 * @param reltol relative tolerance for error control
 * @param abstol absolute tolerance for error control
 * @param hmin minimal step size
 *
 * @return pair of vectors (t_k,y_k)
 */
/* SAM_LISTING_BEGIN_3 */
template <class DiscEvlOp>
std::pair<std::vector<double>, std::vector<double>> OdeIntSsCtrl(
    const DiscEvlOp& Psi, unsigned int p, double y0, double T, double h0,
    double reltol, double abstol, double hmin) {
  double h = h0;
  std::vector<double> Y, t;
  Y.push_back(y0);
  t.push_back(0.0);

  double y = y0;

  while (t.back() < T && h > hmin) {
#if SOLUTION
    double y_high = PsiTilde(Psi, p, std::min(T - t.back(), h), y);
    double y_low = Psi(std::min(T - t.back(), h), y);
    double est = std::abs(y_high - y_low);
    double tol = std::max(reltol * std::abs(y), abstol);

    if (est < tol) {
      y = y_high;
      t.push_back(t.back() + std::min(T - t.back(), h));
      Y.push_back(y_high);
    }
    h = h * std::max(0.5, std::min(2., std::pow(tol / est, 1. / (p + 1))));
#else
    double y_high = 1.0;  // TODO: fix this line
    double y_low = 0.0;   // TODO: fix this line
    double est = std::abs(y_high - y_low);
    if (true /* TODO: fix this line */) {
      y = y_high;
      t.push_back(t.back() + std::min(T - t.back(), h));
      Y.push_back(y_high);
    }
    // TODO: update $h$
#endif
  }
  if (h < hmin) {
    std::cerr << "Warning: Failure at t=" << t.back()
              << ". Unable to meet integration tolerances without reducing the "
                 "step size below the smallest value allowed ("
              << hmin << ") at time t." << std::endl;
  }

  return {t, Y};
}
/* SAM_LISTING_END_3 */

std::pair<std::vector<double>, std::vector<double>> SolveTangentIVP();

}  // namespace ODESolve

#endif
