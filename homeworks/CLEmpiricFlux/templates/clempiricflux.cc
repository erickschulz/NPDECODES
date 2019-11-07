/**
 * @file clempiricflux.cc
 * @brief NPDE exam problem summer 2019 "CLEmpiricFlux" code
 * @author Oliver Rietmann
 * @date 18.07.2019
 * @copyright Developed at ETH Zurich
 */

#include "clempiricflux.h"

#include <cassert>

#include <Eigen/Core>

namespace CLEmpiricFlux {

/**
 * @brief Bisection algorithm for root finding of an increasing function.
 *
 * @param g continuous changing sign in the interval [v, w]
 * @param v lower bound of the interval containg the root
 * @param w upper bound of the interval containg the root
 * @param tol error tolerance for stopping criterion
 * @return approximate root x in [v, w]
 */
template <typename FUNCTOR>
double findRoots(double v, double w, FUNCTOR &&g, double tol = 1.0E-6) {
  double x = v;  // approximate root
  const double len = w - v;
  constexpr static const int maxN = 1000;
  double gv = g(v), gw = g(w);
  // Ensure that function changes sign
  assert(gv * gw <= 0);
  for (int N = 0; (std::abs(w - v) > tol * len) && N < maxN; N++) {
    x = (v + w) / 2.0;
    const double gx = g(x);
    if (gv * gx < 0.0) {
      // Sign change in left half of [v,w]
      w = x;
      gw = gx;
    } else {
      // Sign change in right half of [v,w]
      v = x;
      gv = gx;
    }
  }
  return x;
}

GodunovFlux::GodunovFlux(const UniformCubicSpline &f) : _f(f){};

double GodunovFlux::operator()(double v, double w) const {
  double result;
  /* SOLUTION_BEGIN */
  /* SOLUTION_END */
  return result;
}

}  // namespace CLEmpiricFlux
