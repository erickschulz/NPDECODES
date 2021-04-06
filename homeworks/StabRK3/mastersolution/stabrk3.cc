/**
 * @file stabrk3.cc
 * @brief NPDE homework StabRK3 code
 * @author Unknown, Oliver Rietmann
 * @date 04.04.2021
 * @copyright Developed at ETH Zurich
 */

#include "stabrk3.h"

#include <Eigen/Core>
#include <vector>

#include "rkintegrator.h"

namespace StabRK3 {

/* SAM_LISTING_BEGIN_0 */
Eigen::Vector2d predPrey(Eigen::Vector2d y0, double T, unsigned int N) {
  double h = T / N;
  Eigen::Vector2d y = y0;

  auto f = [](Eigen::Vector2d y) {
    Eigen::Vector2d y_dot;
    y_dot(0) = (1 - y(1)) * y(0);
    y_dot(1) = (y(0) - 1) * y(1);
    return y_dot;
  };

  for (int j = 0; j < N; ++j) {
    Eigen::Vector2d k1 = f(y);
    Eigen::Vector2d k2 = f(y + h * k1);
    Eigen::Vector2d k3 = f(y + (h / 4.) * k1 + (h / 4.) * k2);
    y = y + (h / 6.) * k1 + (h / 6.) * k2 + (2. * h / 3.) * k3;
  }

  return y;
}
/* SAM_LISTING_END_0 */

/* SAM_LISTING_BEGIN_1 */
std::vector<Eigen::Vector2d> simulatePredPrey(
    const std::vector<unsigned int> &N_list) {
  int M = N_list.size();
  std::vector<Eigen::Vector2d> yT_list(M);

  // Initialize RK with Butcher table
  Eigen::Matrix3d A;
  A << 0, 0, 0, 1., 0, 0, 1. / 4., 1. / 4., 0;
  Eigen::Vector3d b(1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0);
  RKIntegrator<Eigen::Vector2d> integrator(A, b);

  // RHS function f for the prey/predator model
  double alpha1 = 1.0;
  double alpha2 = 1.0;
  double beta1 = 1.0;
  double beta2 = 1.0;
  auto f = [&alpha1, &alpha2, &beta1, &beta2](Eigen::Vector2d y) {
    Eigen::Vector2d temp = y;
    temp(0) *= alpha1 - beta1 * y(1);
    temp(1) *= -alpha2 + beta2 * y(0);
    return temp;
  };

  // Compute the solution(s) of the IVP(s)
  double T = 1.0;
  Eigen::Vector2d y0(100.0, 1.0);
  for (int m = 0; m < M; ++m) {
    yT_list[m] = integrator.solve(f, T, y0, N_list[m]).back();
  }

  return yT_list;
}
/* SAM_LISTING_END_1 */

}  // namespace StabRK3
