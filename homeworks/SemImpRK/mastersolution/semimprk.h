#ifndef SEMIMPRK_H_
#define SEMIMPRK_H_

/**
 * @file semimprk.h
 * @brief NPDE homework SemImpRK code
 * @author Unknown, Oliver Rietmann
 * @date 04.04.2021
 * @copyright Developed at ETH Zurich
 */

#include <Eigen/LU>
#include <cmath>
#include <vector>

namespace SemImpRK {

//! \brief Solve the autonomous IVP y' = f(y), y(0) = y0 using Rosenbrock method
//! Use semi-implicit Rosenbrock method using Jacobian evaluation. Equidistant
//! steps of size T/N. \tparam Function function type for r.h.s. f \tparam
//! Jacobian function type for Jacobian df \param[in] f r.h.s. func f \param[in]
//! df Jacobian df of f \param[in] y0 initial data y(0) \param[in] N number of
//! equidistant steps \param[in] T final time \return vector of y_k for each
//! step k from 0 to N
/* SAM_LISTING_BEGIN_0 */
template <class Function, class Jacobian>
std::vector<Eigen::VectorXd> solveRosenbrock(Function &&f, Jacobian &&df,
                                             const Eigen::VectorXd &y0,
                                             unsigned int N, double T) {
  // Will contain all time steps
  std::vector<Eigen::VectorXd> res(N + 1);

  // TO DO: (13-2.c) Implement the Rosenbrock method. Note that the function
  // f should take as argument an Eigen::Eigen::VectorXd, and return one as
  // well. The Jacobian df should take as argument an Eigen::Eigen::VectorXd,
  // but return a square Eigen::Eigen::MatrixXd.
  res.at(0) = y0;  // Push initial data
  const double h = T / N;
  const double a = 1. / (std::sqrt(2) + 2.);

  // Some temporary variables
  Eigen::VectorXd k1, k2;
  Eigen::MatrixXd J, W;

  // Main loop: *up to N (performs N steps)*
  for (unsigned int i = 1; i <= N; ++i) {
    Eigen::VectorXd &yprev = res.at(i - 1);

    // Jacobian computation
    J = df(yprev);
    W = Eigen::MatrixXd::Identity(J.rows(), J.cols()) - a * h * J;

    // Reuse factorization for each step
    auto W_lu = W.partialPivLu();

    // Increments
    k1 = W_lu.solve(f(yprev));
    k2 = W_lu.solve(f(yprev + 0.5 * h * k1) - a * h * J * k1);

    // Push new step
    res.at(i) = yprev + h * k2;
  }
  return res;
}
/* SAM_LISTING_END_0 */

double cvgRosenbrock();

}  // namespace SemImpRK

#endif  // #ifndef SEMIMPRK_H_
