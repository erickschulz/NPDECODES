#ifndef SEMIMPRK_H_
#define SEMIMPRK_H_

/**
 * @file semimprk.h
 * @brief NPDE homework SemImpRK code
 * @author Unknown, Oliver Rietmann
 * @date 04.04.2021
 * @copyright Developed at ETH Zurich
 */

#include <Eigen/Core>
#include <Eigen/LU>
#include <cmath>
#include <vector>

namespace SemImpRK {

/** Solve the autonomous IVP y' =f(y), y(0) = y0 using the Rosenbrock method.*/
/* SAM_LISTING_BEGIN_0 */
template <class Function, class Jacobian>
std::vector<Eigen::VectorXd> SolveRosenbrock(Function &&f, Jacobian &&df,
                                             const Eigen::VectorXd &y0,
                                             unsigned int M, double T) {
  // Will contain all time steps
  std::vector<Eigen::VectorXd> res(M + 1);

  res[0] = y0;  // Push initial data
  const double h = T / M;
  const double a = 1. / (std::sqrt(2) + 2.);

  // Some temporary variables
  Eigen::VectorXd k1, k2;
  Eigen::MatrixXd J, W;

  // Main loop: perform M steps
  for (unsigned int i = 1; i <= M; ++i) {
    Eigen::VectorXd &yprev = res[i - 1];

    // Jacobian computation
    J = df(yprev);
    W = Eigen::MatrixXd::Identity(J.rows(), J.cols()) - a * h * J;

    // Reuse factorization
    auto W_lu = W.partialPivLu();

    // Increments
    k1 = W_lu.solve(f(yprev));
    k2 = W_lu.solve(f(yprev + 0.5 * h * k1) - a * h * J * k1);

    // Push new step
    res[i] = yprev + h * k2;
  }
  return res;
}
/* SAM_LISTING_END_0 */

double CvgRosenbrock();

}  // namespace SemImpRK

#endif  // #ifndef SEMIMPRK_H_
