#ifndef EXPONENTIALINTEGRATOR_H_
#define EXPONENTIALINTEGRATOR_H_

/**
 * @file exponentialintegrator.h
 * @brief NPDE homework ExponentialIntegrator code
 * @author Unknown, Oliver Rietmann
 * @date 04.04.2021
 * @copyright Developed at ETH Zurich
 */

#include <Eigen/Core>

namespace ExponentialIntegrator {

Eigen::MatrixXd phim(const Eigen::MatrixXd &Z);

//! \brief Calculate a single step of the exponential Euler method.
//! \tparam Function function object for r.h.s. function
//! \tparam Jacobian function object for Jacobian of r.h.s.
//! \param[in] y0 The initial state
//! \param[in] f The r.h.s function $f$
//! \param[in] df The Jacobian of $f$
//! \param[in] h The stepsize of the method.
//! \return A single step of the Exponential Euler method
/* SAM_LISTING_BEGIN_0 */
template <class Function, class Jacobian>
Eigen::VectorXd exponentialEulerStep(const Eigen::VectorXd &y0, Function &&f,
                                     Jacobian &&df, double h) {
#if SOLUTION
  return y0 + h * phim(h * df(y0)) * f(y0);
#else
  //====================
  // Your code goes here
  // Replace the following dummy return value:
  return y0;
  //====================
#endif
}
/* SAM_LISTING_END_0 */

void testExpEulerLogODE();

}  // namespace ExponentialIntegrator

#endif  // #ifndef EXPONENTIALINTEGRATOR_H_
