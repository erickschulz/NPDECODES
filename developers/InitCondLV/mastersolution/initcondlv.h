#ifndef LV_H_
#define LV_H_
////
//// Copyright (C) 2016 SAM (D-MATH) @ ETH Zurich
//// Author(s): lfilippo <filippo.leonardi@sam.math.ethz.ch>
//// Contributors: tille, jgacon, dcasati
//// This file is part of the NumCSE repository.

#include <Eigen/Core>
#include <iostream>
#include <utility>

#include "ode45.h"

namespace InitCondLV {

/*!
 * \brief Compute the maps Phi and W at time T.
 * Use initial data given by u0 and v0.
 * \param[in] u0 First component.
 * \param[in] v0 Second component.
 * \param[in] T Final time.
 */
/* SAM_LISTING_BEGIN_1 */
std::pair<Eigen::Vector2d, Eigen::Matrix2d> PhiAndW(double u0, double v0,
                                                    double T) {
  std::pair<Eigen::Vector2d, Eigen::Matrix2d> PaW;
  // TO DO: (12-7.f)  Calculate the time evolution of Phi(t,y0) and W(t,y0)
  // up to the final time T, using the inital value y0=[u0,v0].
  // Save the values of Phi and W at time T in PaW.first and
  // PaW.second respectively.
  // Hint: First write the ODE
  // [ y'(t), W'(t,y0) ] = [ f(y(t)), Df(y(t))W(t,y0) ]
  // in vectorized form, i.e. find a function g:R^6 -> R^6
  // such that the ODE w'=g(w) is equivalent.
  // Then use the ode45 class for the function g.
#if SOLUTION
  auto f = [](const Eigen::VectorXd& w) {
    Eigen::VectorXd temp(6);
    temp(0) = (2. - w(1)) * w(0);
    temp(1) = (w(0) - 1.) * w(1);
    temp(2) = (2. - w(1)) * w(2) - w(0) * w(3);
    temp(3) = w(1) * w(2) + (w(0) - 1.) * w(3);
    temp(4) = (2. - w(1)) * w(4) - w(0) * w(5);
    temp(5) = w(1) * w(4) + (w(0) - 1.) * w(5);
    return temp;
  };

  Eigen::VectorXd w0(6);
  w0 << u0, v0, 1., 0, 0, 1.;

  // Construct ode solver with r.h.s
  ode45<Eigen::VectorXd> O(f);
  // Set options
  O.options.rtol = 1e-14;
  O.options.atol = 1e-12;
  // Solve ODE
  auto sol = O.solve(w0, T);
  // Extract needed component
  Eigen::VectorXd wT = sol.back().first;

  PaW.first << wT(0), wT(1);
  PaW.second << wT(2), wT(4), wT(3), wT(5);
#else
  //====================
  // Your code goes here
  //====================
#endif
  return PaW;
}
/* SAM_LISTING_END_1 */

/*!
 * \brief Compute initial conditions u0, v0 such that solution has period T=5
 * inital approximation [3,2]^T
 * \param[out] u0
 * \param[out] v0
 */
/* SAM_LISTING_BEGIN_2 */
Eigen::Vector2d findInitCond(void) {
  Eigen::Vector2d y;
  // TO DO: (12-7.g) Use Newton's method with initial guess (3,2) to find a
  // zero of the function F from (12-7.c) for the period $T_P=5$.
  // Store the resulting vector in y.
  // Hint: F and its Jacobian matrix can be computed from the result of
  // PhiAndW().
#if SOLUTION
  y << 3, 2;     // Initial guess
  double T = 5;  // Period we require

  // Compute Phi and W from initial guess
  std::pair<Eigen::Vector2d, Eigen::Matrix2d> PaW = PhiAndW(y(0), y(1), T);
  Eigen::Vector2d F = PaW.first - y;  // Value of F at initial guess
  Eigen::Matrix2d DF;                 // Declare Jacobian of F

  // Until we are happy with our approximation
  while (F.norm() > 1e-5) {  // Residual based termination
    // Calculate Jacobian
    DF = PaW.second - Eigen::MatrixXd::Identity(2, 2);
    // Use Newton iteration
    y = y - DF.lu().solve(F);
    // Test current guess
    PaW = PhiAndW(y(0), y(1), T);
    // Get value of F at the current guess
    F = PaW.first - y;
  }
#else
  //====================
  // Your code goes here
  //====================
#endif

  // TO DO: (12-7.g) Calculate the evolution of the Lotka-Volterra ODE up to
  // time 100, using the initial value y, found by the Newton iterations above.
  // If the evolution is indeed 5-periodic, then the solution at time 100
  // should have the value y. Print a warning message if this is not the case.
#if SOLUTION
  PaW = PhiAndW(y(0), y(1), 100);
  if ((y - PaW.first).norm() > 1e-5) {
    std::cout << "Warning: Solution not periodic, y(100) != y(0)" << std::endl;
  }
#else
  //====================
  // Your code goes here
  //====================
#endif

  return y;
}
/* SAM_LISTING_END_2 */

}  // namespace InitCondLV

#endif  // #define LV_H_
