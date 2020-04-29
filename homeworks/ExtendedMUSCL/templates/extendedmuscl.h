#ifndef EXTENDEDMUSCL_H_
#define EXTENDEDMUSCL_H_

/**
 * @file extendedmuscl.h
 * @brief NPDE homework ExtendedMUSCL code
 * @author Oliver Rietmann
 * @date 04.08.2019
 * @copyright Developed at ETH Zurich
 */

#include "slopelimfluxdiff.h"

#include <cassert>
#include <cmath>

#include <Eigen/Core>

namespace ExtendedMUSCL {

/**
 * @brief Computes the Godunov numerical Flux for
 * flux function f(u) = u(log(u) - 1) and u > 0.
 *
 * @param v left state, v > 0
 * @param w right state, w > 0
 * @return F_GD(v, w).
 */
double logGodunovFlux(double v, double w);

/**
 * @brief Computes the scaled slopes for linear reconstruction,
 * so that it serves as argument for slopelimfluxdiff(per).
 *
 * @param mu_left $\mu_{j-1}$
 * @param mu_center $\mu_{j}$
 * @param mu_right $\mu_{j+1}$
 * @return MC slope limiter as in the problem description
 */
double limiterMC(double mu_left, double mu_center, double mu_right);

/**
 * @brief Solves the ODE $\dot{y} = f(y)$ using the SSP
 * method given in the problem description.
 *
 * @param f right-hand side of ODE modeling std::function<State(State)>
 * @param y inital data, State supports addition and multiplication
 * by a double (from the left). Examples: Eigen::VectorXd or double
 * @param tau small timestep tau > 0
 * @return solution y after time tau.
 */
/* SAM_LISTING_BEGIN_1 */
template <typename FUNCTOR, typename State>
State sspEvolop(FUNCTOR &&f, State y, double tau) {
  State y_tau;
  //====================
  // Your code goes here
  //====================
  return y_tau;
}
/* SAM_LISTING_END_1 */

/**
 * @brief Solves $\dot{u} + \partial_xf(u) = 0$ for
 * $f(u)=u(\log(u)-1)$ using the finite volume method
 * from the problem description.
 *
 * @param u0 inital data modeling std::function<double(double)>
 * @param T final time, T > 0
 * @param n number of finite volume cells
 * @return approximate solution u(x, T).
 */
/* SAM_LISTING_BEGIN_4 */
template <typename U0_FUNCTOR>
Eigen::VectorXd solveClaw(U0_FUNCTOR &&u0, double T, unsigned int n) {
  // Set up spacial mesh and inital data.
  double a = 0.0;
  double b = 1.0;
  double h = (b - a) / n;
  Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(n, a + 0.5 * h, b - 0.5 * h);

  // Approximate dual cell averages at t=0
  Eigen::VectorXd mu = x.unaryExpr(u0);

  double alpha = mu.minCoeff(); // lower bound for initial data
  double beta = mu.maxCoeff();  // upper bound for initial data
  assert(alpha > 0.0 && beta > 0.0);

  //====================
  // Your code goes here
  //====================

  return mu;
}
/* SAM_LISTING_END_4 */

} // namespace ExtendedMUSCL

#endif // EXTENDEDMUSCL_H_
