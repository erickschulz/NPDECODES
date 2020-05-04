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
KineticPropagator::KineticPropagator(const SparseMatrixXd &A, const SparseMatrixXcd &M, double tau) {
  B_plus_ = M + 0.5 * tau * A.cast<std::complex<double>>();
  SparseMatrixXcd B_minus = M - 0.5 * tau * A.cast<std::complex<double>>();
  solver_.compute(B_minus);
}

Eigen::VectorXcd KineticPropagator::operator()(const Eigen::VectorXcd &mu) const {
  return solver_.solve(B_plus_ * mu);
}

// InteractionPropagator
InteractionPropagator::InteractionPropagator(double tau) {
  phase_multiplier_ = [tau] (std::complex<double> z) {
    const std::complex<double> i(0, 1);
    return std::exp(-i * tau * std::norm(z)) * z;
  };
}

Eigen::VectorXcd InteractionPropagator::operator()(const Eigen::VectorXcd &mu) const {
  return mu.unaryExpr(phase_multiplier_);
}

}  // namespace NonLinSchroedingerEquation
