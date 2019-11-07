#ifndef UNIFORMCUBICSPLINE_H
#define UNIFORMCUBICSPLINE_H

/**
 * @file uniformcubicspline.h
 * @brief NPDE exam problem summer 2019 "CLEmpiricFlux" code
 * @author Oliver Rietmann
 * @date 18.07.2019
 * @copyright Developed at ETH Zurich
 */

#include <Eigen/Core>

namespace CLEmpiricFlux {

class UniformCubicSpline {
 public:
  UniformCubicSpline(double a, double b, const Eigen::VectorXd f,
                     const Eigen::VectorXd M);
  double operator()(double u) const;  // Point evaluation operator
  double derivative(double u) const;  // Evaluation of derivative
 private:
  unsigned int _n;     // Number of nodes - 1
  double _a, _b;       // Interval boundaries
  Eigen::VectorXd _f;  // Values of flux function at nodes
  Eigen::VectorXd _M;  // Values of second derivatives of f at nodes
};

}  // namespace CLEmpiricFlux

#endif
