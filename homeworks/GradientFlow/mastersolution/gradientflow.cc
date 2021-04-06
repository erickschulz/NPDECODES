/**
 * @file gradientflow.h
 * @brief NPDE homework GradientFlow code
 * @author Unknown, Oliver Rietmann
 * @date 04.04.2021
 * @copyright Developed at ETH Zurich
 */

#include "gradientflow.h"

#include <Eigen/Core>
#include <cmath>
#include <vector>

namespace GradientFlow {

/* SAM_LISTING_BEGIN_0 */
Eigen::MatrixXd ButcherMatrix() {
  Eigen::MatrixXd A(6, 5);
  // clang-format off
  A <<       0.25,          0.,       0.,       0.,   0.,
    0.5,        0.25,       0.,       0.,   0.,
    17./50.,     -1./25.,     0.25,       0.,   0.,
    371./1360., -137./2720., 15./544.,     0.25,   0.,
    25./24.,    -49./48., 125./16., -85./12., 0.25,
    25./24.,    -49./48., 125./16., -85./12., 0.25;
  // clang-format on
  return A;
}
/* SAM_LISTING_END_0 */

/* SAM_LISTING_BEGIN_1 */
std::vector<Eigen::VectorXd> solveGradientFlow(const Eigen::VectorXd &d,
                                               double lambda,
                                               const Eigen::VectorXd &y,
                                               double T, unsigned int N) {
  std::vector<Eigen::VectorXd> Y(N + 1);
  // TO DO (0-2.h)
  // Define the right hand side of the ODE y' = f(y), and the Jacobian of f.
  auto f = [d, lambda](const Eigen::VectorXd &yt) {
    Eigen::VectorXd val =
        -2. * std::cos(yt.squaredNorm()) * yt - 2. * lambda * d.dot(yt) * d;
    return val;
  };
  auto J = [d, lambda](const Eigen::VectorXd &yt) {
    int dim = yt.size();
    Eigen::MatrixXd term1 =
        4. * std::sin(yt.squaredNorm()) * yt * yt.transpose();
    Eigen::MatrixXd term2 =
        -2. * std::cos(yt.squaredNorm()) * Eigen::MatrixXd::Identity(dim, dim);
    Eigen::MatrixXd term3 = -2. * lambda * d * d.transpose();
    Eigen::MatrixXd Jval = term1 + term2 + term3;
    return Jval;
  };

  // Split the interval [0,T] into N intervals of size h.
  double h = T / N;
  Eigen::VectorXd yt = y;
  Y.at(0) = y;
  // Evolve up to time T:
  for (int i = 1; i <= N; i++) {
    yt = discEvolSDIRK(f, J, yt, h);
    Y.at(i) = yt;
  }
  return Y;
}
/* SAM_LISTING_END_1 */

}  // namespace GradientFlow
