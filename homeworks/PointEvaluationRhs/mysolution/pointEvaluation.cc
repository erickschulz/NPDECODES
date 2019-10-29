/**
 * @ file pointEvaluation.cc
 * @ brief NPDE homework PointEvaluationRhs code
 * @ author Christian Mitsch
 * @ date 22.03.2019
 * @ copyright Developed at ETH Zurich
 */

#include "pointEvaluation.h"
#include <cmath>
#include <iostream>
#include "norms.h"

namespace PointEvaluationRhs {

/* SAM_LISTING_BEGIN_5 */
Eigen::Vector2d GlobalInverseTria(Eigen::Matrix<double, 2, 3> mycorners,
                                  Eigen::Vector2d x) {
  Eigen::Vector2d x_hat;

  /* BEGIN_SOLUTION */
  /* TODO Your implementation goes here! */
  /* END_SOLUTION */
  return x_hat;
}
/* SAM_LISTING_END_5 */

/** @brief Numerically stable solution of a quadratic equation in R
 * @param a,b,c coefficients of quadratic polynomias ax^2+bx+c
 * @return both zeros, NaN if complex
 */
/* SAM_LISTING_BEGIN_4 */
std::pair<double, double> solveQuadraticEquation(double a, double b, double c) {
  // Implement the cases which are solvable and return their solutions
  /* BEGIN_SOLUTION */
  /* TODO Your implementation goes here! */
  /* END_SOLUTION */
  // Return NAN if there are no (real) roots
  return {NAN, NAN};
}

/* SAM_LISTING_END_4 */

/** @brief Computes the area of a triangle
 * @param a,b,c vertex coordinate vectors
 */
inline double triaArea(const Eigen::Vector2d a, const Eigen::Vector2d b,
                       const Eigen::Vector2d c) {
  double result = 0;
  /* BEGIN_SOLUTION */
  /* TODO Your implementation goes here! */
  /* END_SOLUTION */
  return result;
}

Eigen::Vector2d GlobalInverseQuad(Eigen::Matrix<double, 2, 4> vert,
                                  Eigen::Vector2d x) {
  constexpr double kEPS = 1.0E-8;
  Eigen::Vector2d x_hat;

  // Implement and use the functions triaArea and solveQuadraticEquation
  /* BEGIN_SOLUTION */
  /* TODO Your implementation goes here! */
  /* END_SOLUTION */
  return x_hat;
}

std::pair<double, double> normsSolutionPointLoadDirichletBVP(
      const lf::assemble::DofHandler &dofh, Eigen::Vector2d source_point,
      Eigen::VectorXd &sol_vec) {
  std::pair<double, double> result(0, 0);
  const unsigned int N_dofs = dofh.NoDofs();
  sol_vec.resize(N_dofs);
  sol_vec.setZero();
  /* BEGIN_SOLUTION */
  /* TODO Your implementation goes here! */
  /* END_SOLUTION */
  return result;
}

/* SAM_LISTING_BEGIN_6 */
Eigen::VectorXd DeltaLocalVectorAssembler::Eval(const lf::mesh::Entity &cell) {
  Eigen::VectorXd result;
  // get the coordinates of the corners of this cell
  lf::geometry::Geometry *geo_ptr = cell.Geometry();
  auto vertices = lf::geometry::Corners(*geo_ptr);
  /* BEGIN_SOLUTION */
  /* TODO Your implementation goes here! */
  /* END_SOLUTION */
  return result;
}
/* SAM_LISTING_END_6 */

}  // namespace PointEvaluationRhs
