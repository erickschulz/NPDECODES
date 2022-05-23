#ifndef SEMILAGRANGIAN_H
#define SEMILAGRANGIAN_H
/**
 * @file semilagrangian.h
 * @brief NPDE homework TEMPLATE HEADER FILE
 * @author
 * @date May 2022
 * @copyright Developed at SAM, ETH Zurich
 */

#include <Eigen/Core>
#include <cmath>
#include <iostream>

namespace SemiLagrangian {

/* SAM_LISTING_BEGIN_1 */
/**
 * @brief Evaluates the solution of the transport problem using the solution
 * formula and Heun's method
 * @param x location
 * @param K maximum number of steps of Heun's method
 * @param t time
 * @param v velocity field (time independent)
 * @param u0 initial condition
 * @return u(x,t)
 */
template <typename FUNCTOR_V, typename FUNCTOR_U0>
double solveTransport(const Eigen::Vector2d& x, int K, double t, FUNCTOR_V&& v,
                      FUNCTOR_U0&& u0) {
#if SOLUTION
  double tau = t / K;     // timestep
  Eigen::Vector2d y = x;  // starting point

  for (int i = 0; i < K; ++i) {
    // A single step of Heun's method, expressed by RK increments
    Eigen::Vector2d k_1 = v(y);
    Eigen::Vector2d k_2 = v(y - 2. / 3. * tau * k_1);
    y -= tau / 4. * k_1 + 3. / 4. * tau * k_2;

    // check, if current point is outside the domain
    if (y(0) < 0. || y(0) > 1. || y(1) < 0. || y(1) > 1.) {
      return 0.0;
    }
  }
  return u0(y);
#else
  //====================
  // Your code goes here
  //====================
  return 0.0;
#endif
}
/* SAM_LISTING_END_1 */
/**
 * @brief Evaluates an FE function u_h
 * @param x point coordinates
 * @param u FE basis coefficients of u_h
 * @return u_h(x)
 */
double evalFEfunction(const Eigen::Vector2d& x, const Eigen::VectorXd& u);

/**
 * @brief Computes the interior nodes of the mesh
 * @param M number of cells in one direction
 * @return 2x(M-1)^2 matrix containing the interior nodes of the mesh as columns
 */
Eigen::MatrixXd findGrid(int M);

/* SAM_LISTING_BEGIN_2 */
/**
 * @brief Evaluates the source term of the semi-lagrangian scheme.
 * @param u_old FE basis coefficients of the solution at the previous timestep
 * @param tau step size
 * @param velocity velociy field
 * @return RHS-vector of the scheme
 */
template <typename FUNCTOR>
Eigen::VectorXd semiLagrangeSource(const Eigen::VectorXd& u_old, double tau,
                                   FUNCTOR&& velocity) {
  // Note: components of coefficient vectors are associated
  // with interior nodes only
  int N = u_old.size();
  // Extract number of cells in one direction
  int root = std::round(std::sqrt(N));
  if (N != root * root) {
    std::cerr << "The number of dofs should be a perfect square!" << std::endl;
  }
  int M = root + 1;
  Eigen::MatrixXd grid = findGrid(M);
  Eigen::VectorXd f(N);

#if SOLUTION
  for (int i = 0; i < grid.cols(); ++i) {
    // Find grid point corresponding to a degree of freedom
    Eigen::Vector2d x = grid.col(i);
    // Determine location of advected gridpoint
    Eigen::Vector2d y = x - tau * velocity(x);
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
  // Finally scale with $h^{2}$
  return f / (M * M);  // * 1 (from $[0,1]^2$) * 4 (from no. of adjacent
                       // squares) / 4 (from no. of vertices of square)
#else
  //====================
  // Your code goes here
  //====================
  return f;
#endif
}
/* SAM_LISTING_END_2 */

/**
 * @brief Solves the pure transport problem using the semi-lagrangian scheme.
 * @param M number of cells in one direction
 * @param K total number of equidistant timesteps
 * @param T final time
 * @return FE Basis coefficients of the approximate solution at final time
 */
Eigen::VectorXd semiLagrangePureTransport(int M, int K, double T);

void testFloorAndDivision();

}  // namespace SemiLagrangian

#endif  // SEMILAGRANGIAN_H