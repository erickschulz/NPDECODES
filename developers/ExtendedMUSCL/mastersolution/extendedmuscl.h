#ifndef EXTENDEDMUSCL_H_
#define EXTENDEDMUSCL_H_

/**
 * @file extendedmuscl.h
 * @brief NPDE homework ExtendedMUSCL code
 * @author Oliver Rietmann
 * @date 04.08.2019
 * @copyright Developed at ETH Zurich
 */

#include <Eigen/Core>
#include <cassert>
#include <cmath>

#include "slopelimfluxdiff.h"

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
#if SOLUTION
  State k = y + tau * f(y);
  k = 0.75 * y + 0.25 * k + 0.25 * tau * f(k);
  y_tau = 1.0 / 3.0 * y + 2.0 / 3.0 * k + 2.0 / 3.0 * tau * f(k);
#else
  //====================
  // Your code goes here
  //====================
#endif
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

  double alpha = mu.minCoeff();  // lower bound for initial data
  double beta = mu.maxCoeff();   // upper bound for initial data
  assert(alpha > 0.0 && beta > 0.0);

#if SOLUTION
  // Set timestep tau according to the CFL-condition
  double tau =
      h / std::max(std::abs(std::log(alpha)), std::abs(std::log(beta)));
  // Semi-discretization (discretization in space)
  double h_inv = 1.0 / h;
  // Define right-hand-side for the SSP ODE solver
  auto semi_discrete_rhs =
      [h_inv](const Eigen::VectorXd &mu) -> Eigen::VectorXd {
    return -h_inv * slopelimfluxdiffper(mu, &logGodunovFlux, &limiterMC);
  };
  // Timestepping: Solve the semi-discrete ODE
  int N = (int)(T / tau + 0.5);
  for (int i = 0; i < N; ++i) mu = sspEvolop(semi_discrete_rhs, mu, tau);
#else
  //====================
  // Your code goes here
  //====================
#endif

  return mu;
}
/* SAM_LISTING_END_4 */

}  // namespace ExtendedMUSCL

#endif  // EXTENDEDMUSCL_H_
