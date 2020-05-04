/**
 * @file propagator.h
 * @brief NPDE homework NonLinSchroedingerEquation code
 * @author Oliver Rietmann
 * @date 04.05.2020
 * @copyright Developed at ETH Zurich
 */

#include "propagator.h"

#include <cmath>
#include <complex>

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>

namespace NonLinSchroedingerEquation {

// KineticPropagator
/* SAM_LISTING_BEGIN_1 */
KineticPropagator::KineticPropagator(const SparseMatrixXd &A, const SparseMatrixXcd &M, double tau) {
#if SOLUTION
  B_plus_ = M + 0.5 * tau * A.cast<std::complex<double>>();
  SparseMatrixXcd B_minus = M - 0.5 * tau * A.cast<std::complex<double>>();
  solver_.compute(B_minus);
#else
  //====================
  // Your code goes here
  //====================
#endif
}

Eigen::VectorXcd KineticPropagator::operator()(const Eigen::VectorXcd &mu) const {
#if SOLUTION
  return solver_.solve(B_plus_ * mu);
#else
  //====================
  // Your code goes here
  // Replace mu by its value after a timestep tau
  return mu;
  //====================
#endif
}
/* SAM_LISTING_END_1 */

// InteractionPropagator
/* SAM_LISTING_BEGIN_2 */
InteractionPropagator::InteractionPropagator(double tau) {
#if SOLUTION
  phase_multiplier_ = [tau] (std::complex<double> z) {
    const std::complex<double> i(0, 1);
    return std::exp(-i * tau * std::norm(z)) * z;
  };
#else
  //====================
  // Your code goes here
  // Replace this by the proper phase multiplier
  return [] (std::complex<double> z) { return z; };
  //====================
#endif
}

Eigen::VectorXcd InteractionPropagator::operator()(const Eigen::VectorXcd &mu) const {
#if SOLUTION
  return mu.unaryExpr(phase_multiplier_);
#else
  //====================
  // Your code goes here
  // Replace mu by its value after a timestep tau
  return mu;
  //====================
#endif
}
/* SAM_LISTING_END_2 */

}  // namespace NonLinSchroedingerEquation
