#ifndef PROPAGATOR_H_
#define PROPAGATOR_H_

/**
 * @file propagator.h
 * @brief NPDE homework NonLinSchroedingerEquation code
 * @author Oliver Rietmann
 * @date 03.05.2020
 * @copyright Developed at ETH Zurich
 */

#include <complex>

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>

namespace NonLinSchroedingerEquation {

class Propagator {
public:
  Propagator() = default;
  virtual ~Propagator() = default;
  virtual Eigen::VectorXcd operator()(const Eigen::VectorXcd &mu) const = 0;

private:
  Propagator(const Propagator &) = delete;
  Propagator(Propagator &&) = delete;
  Propagator & operator=(const Propagator &) = delete;
  Propagator & operator =(const Propagator &&) = delete;
};

class KineticPropagator : public Propagator {
  using SparseMatrixXcd = Eigen::SparseMatrix<std::complex<double>>;
  using SparseMatrixXd = Eigen::SparseMatrix<double>;

public:
  KineticPropagator(const SparseMatrixXd &A, const SparseMatrixXcd &M, double tau);
  Eigen::VectorXcd operator()(const Eigen::VectorXcd &mu) const override;

private:
  SparseMatrixXcd B_plus_;
  Eigen::SparseLU<SparseMatrixXcd> solver_;
};

class InteractionPropagator : public Propagator {
public:
  InteractionPropagator(double tau);
  Eigen::VectorXcd operator()(const Eigen::VectorXcd &mu) const override;

private:
  std::function<std::complex<double>(std::complex<double>)> phase_multiplier_;
};

}  // namespace NonLinSchroedingerEquation

#endif  // PROPAGATOR_H_
