/**
 * @file finitevolumesineconslaw.cc
 * @brief NPDE homework "FiniteVolumeSineConsLaw" code
 * @author Oliver Rietmann
 * @date 25.04.2019
 * @copyright Developed at ETH Zurich
 */

#include "finitevolumesineconslaw.h"

#include <cmath>

#include <Eigen/Core>

namespace FiniteVolumeSineConsLaw {

constexpr double PI = 3.14159265358979323846;

double f(double x) { return std::sin(PI * x); }

double sineGodFlux(double v, double w) {
  double result;
  /* BEGIN_SOLUTION */
  /* TODO Your implementation goes here! */
  /* END_SOLUTION */
  return result;
}

Eigen::VectorXd sineClawRhs(const Eigen::VectorXd &mu) {
  int N = mu.size();
  double h = 12.0 / N;
  Eigen::VectorXd result(N);

  /* BEGIN_SOLUTION */
  /* TODO Your implementation goes here! */
  /* END_SOLUTION */

  return result;
}

bool blowup(const Eigen::VectorXd &mu) {
  return mu.minCoeff() < 0.0 || mu.maxCoeff() > 2.0;
}

unsigned int findTimesteps() {
  const unsigned int N = 600;
  unsigned int ML = 100;
  unsigned int MR = 200;

  /* BEGIN_SOLUTION */
  /* TODO Your implementation goes here! */
  /* END_SOLUTION */

  return MR;
}

Eigen::VectorXd sineClawReactionRhs(const Eigen::VectorXd &mu, double c) {
  Eigen::VectorXd rhs(mu.size());
  /* BEGIN_SOLUTION */
  /* TODO Your implementation goes here! */
  /* END_SOLUTION */
  return rhs;
}

}  // namespace FiniteVolumeSineConsLaw
