#ifndef SEMILAGRANGIAN_H
#define SEMILAGRANGIAN_H
/**
 * @ file semilagrangian.h
 * @ brief NPDE homework TEMPLATE HEADER FILE
 * @ author
 * @ date May 2022
 * @ copyright Developed at SAM, ETH Zurich
 */

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <cmath>
#include <iostream>
#include <limits>


namespace SemiLagrangian {

/* SAM_LISTING_BEGIN_1 */
template <typename FUNCTOR_V, typename FUNCTOR_U0>
double solveTransport(const Eigen::Vector2d& x, int K, double t, FUNCTOR_V&& v,
                      FUNCTOR_U0&& u0) {
  double tau = t / K;  // timestep
  Eigen::Vector2d y = x;       // starting point

  for (int i = 0; i < K; ++i) {
    // A single step of Heun's method, expressed by RK increments
    Eigen::Vector2d k_1 = v(y);
    Eigen::Vector2d k_2 = v(y - 2. / 3. * tau * k_1);
    y -= tau / 4. * k_1 + 3. / 4. * tau * k_2;
  }

  // Test whether trajectory has hit inflow boundary
  if (y(0) >= 0. && y(0) <= 1. && y(1) >= 0. && y(1) <= 1.) {
    return u0(y);
  } else {
    return 0.;
  }
}
/* SAM_LISTING_END_1 */

double evalFEfunction(const Eigen::Vector2d& x, const Eigen::VectorXd& u);

Eigen::MatrixXd findGrid(int M);

/* SAM_LISTING_BEGIN_2 */
template <typename FUNCTOR>
Eigen::VectorXd semiLagrangeSource(const Eigen::VectorXd& u_old, double tau, double t,
                            FUNCTOR&& velocity) {
  // Note: components of coefficient vectors are associated
  // with interior nodes only
  int N = u_old.size();  // assume dofs on boundary already removed
  // Extract number of cells in one direction
  int root = std::round(std::sqrt(N));
  if (N != root * root) {
    std::cerr << "The number of dofs should be a perfect square!" << std::endl;
  }
  int M = root + 1;

  Eigen::MatrixXd grid = findGrid(M);
  Eigen::VectorXd f(N);
  for (int i = 0; i < grid.cols(); ++i) {
    // Find grid point corresponding to a degree of freedom
    Eigen::Vector2d x = grid.col(i);
    // Determine location of advected gridpoint
    Eigen::Vector2d y = x - tau * velocity(t, x);
    // Test whether it still lies inside the computational domain
    if (y(0) >= 0. && y(0) <= 1. && y(1) >= 0. && y(1) <= 1.) {
      // Evaluate finite element function from previous timestep
      // at preimage of gridpoint under flow
      f(i) = evalFEfunction(y, u_old);
    } else {
      // Zero, if advected point outside domain
      f(i) = 0.;
    }
  }
  // Finally scale with $h^{-2}$
  return f / (M * M);  // * 1 (from $[0,1]^2$) * 4 (from no. of adjacent
                       // squares) / 4 (from no. of vertices of square)
}
/* SAM_LISTING_END_2 */


Eigen::VectorXd semiLagrangePureTransport(int M, int K, double T);

} //namespace SemiLagrangian

#endif //SEMILAGRANGIAN_H