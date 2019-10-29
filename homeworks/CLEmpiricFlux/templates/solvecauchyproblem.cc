/**
 * @file solvecauchyproblem.cc
 * @brief NPDE exam problem summer 2019 "CLEmpiricFlux" code
 * @author Oliver Rietmann
 * @date 19.07.2019
 * @copyright Developed at ETH Zurich
 */

#include "solvecauchyproblem.h"
#include "uniformcubicspline.h"

#include <cmath>

#include <Eigen/Core>

namespace CLEmpiricFlux {

Eigen::Vector2d findSupport(const UniformCubicSpline &f,
                            const Eigen::Vector2d &initsupp, double t) {
  Eigen::Vector2d result;
  /* SOLUTION_BEGIN */
  /* SOLUTION_END */
  return result;
}

template <typename FUNCTOR>
Eigen::VectorXd RalstonODESolver(FUNCTOR &&rhs, Eigen::VectorXd mu0, double tau,
                                 int n) {
  /* SOLUTION_BEGIN */
  /* SOLUTION_END */
  return mu0;
}

Eigen::VectorXd solveCauchyProblem(const UniformCubicSpline &f,
                                   const Eigen::VectorXd &mu0, double h,
                                   double T) {
  Eigen::VectorXd muT(mu0.size());
  /* SOLUTION_BEGIN */
  /* SOLUTION_END */
  return muT;
}

} // namespace CLEmpiricFlux
