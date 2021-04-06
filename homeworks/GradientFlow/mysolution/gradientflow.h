#ifndef GRADIENTFLOW_H_
#define GRADIENTFLOW_H_

/**
 * @file gradientflow.h
 * @brief NPDE homework GradientFlow code
 * @author Unknown, Oliver Rietmann
 * @date 04.04.2021
 * @copyright Developed at ETH Zurich
 */

#include <Eigen/LU>
#include <array>
#include <vector>

namespace GradientFlow {

Eigen::MatrixXd ButcherMatrix();

/* SAM_LISTING_BEGIN_0 */
template <typename Functor, typename Jacobian>
Eigen::VectorXd solveGenStageEquation(Functor &&f, Jacobian &&J,
                                      const Eigen::VectorXd &y,
                                      const Eigen::VectorXd &b, double h,
                                      double rtol = 1E-6, double atol = 1E-8) {
  // Need to solve the equation lhs(g) = g - h*f(y+g)/4 - b = 0.
  auto lhs = [f, y, b, h](const Eigen::VectorXd &g) {
    // lhs(g) = g - h*f(y+g)/4 - b
    Eigen::VectorXd val = g - 0.25 * h * f(y + g) - b;
    return val;
  };
  auto Jlhs = [J, y, h](const Eigen::VectorXd &g) {
    // lhs(g) = g - h*f(y+g)/4 - b, so the Jacobian is
    // Jlhs(g) = Id - h*Jf(y+g)/4
    int dim = y.size();
    Eigen::MatrixXd Jval =
        Eigen::MatrixXd::Identity(dim, dim) - 0.25 * h * J(y + g);
    return Jval;
  };

  // Perform Newton iterations:
  Eigen::VectorXd g = Eigen::VectorXd::Zero(y.size());  // initial guess g=0.
  Eigen::VectorXd delta =
      -Jlhs(g).lu().solve(lhs(g));  // Newton correction term.
  int iter = 0,
      maxiter = 100;  // If correction based termination does not work.
  while (delta.norm() > atol && delta.norm() > rtol * g.norm() &&
         iter < maxiter) {
    g = g + delta;
    delta = -Jlhs(g).lu().solve(lhs(g));
    iter++;
  }
  return g + delta;  // Perform the final step
}
/* SAM_LISTING_END_0 */

/* SAM_LISTING_BEGIN_1 */
template <typename Functor, typename Jacobian>
std::array<Eigen::VectorXd, 5> computeStages(Functor &&f, Jacobian &&J,
                                             const Eigen::VectorXd &y, double h,
                                             double rtol = 1E-6,
                                             double atol = 1E-8) {
  std::array<Eigen::VectorXd, 5> G;  // array of stages
  int d = y.size();
  Eigen::MatrixXd Coeffs = ButcherMatrix();
  // TO DO (0-2.d)
  //====================
  // Your code goes here
  //====================
  return G;
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
template <typename Functor, typename Jacobian>
Eigen::VectorXd solveGenIncrementEquation(Functor &&f, Jacobian &&J,
                                          const Eigen::VectorXd &y,
                                          const Eigen::VectorXd &b, double h,
                                          double rtol = 1E-6,
                                          double atol = 1E-8) {
  // Need to solve the equation lhs(k) = k - f(y + h*k/4 + b) = 0
  auto lhs = [f, y, b, h](const Eigen::VectorXd &k) {
    // lhs(k) = k - f(y + h*k/4 + b)
    Eigen::VectorXd val = k - f(y + 0.25 * h * k + b);
    return val;
  };
  auto Jlhs = [J, y, b, h](const Eigen::VectorXd &k) {
    // lhs(k) = k - f(y + h*k/4 + b), so the Jacobian is
    // Jlhs(g) = Id - Jf(y+h*k/4+b)*h/4
    int dim = y.size();
    Eigen::MatrixXd Jval = Eigen::MatrixXd::Identity(dim, dim) -
                           0.25 * h * J(y + 0.25 * h * k + b);
    return Jval;
  };
  // Perform Newton iterations:
  Eigen::VectorXd k = Eigen::VectorXd::Zero(y.size());  // initial guess k=0.
  Eigen::VectorXd delta =
      -Jlhs(k).lu().solve(lhs(k));  // Newton correction term.
  int iter = 0,
      maxiter = 100;  // If correction based termination does not work.
  while (delta.norm() > atol && delta.norm() > rtol * k.norm() &&
         iter < maxiter) {
    k = k + delta;
    delta = -Jlhs(k).lu().solve(lhs(k));
    iter++;
  }
  return k + delta;  // Perform the final step
}
/* SAM_LISTING_END_2 */

/* SAM_LISTING_BEGIN_3 */
template <typename Functor, typename Jacobian>
std::array<Eigen::VectorXd, 5> computeIncrements(Functor &&f, Jacobian &&J,
                                                 const Eigen::VectorXd &y,
                                                 double h, double rtol = 1E-6,
                                                 double atol = 1E-8) {
  std::array<Eigen::VectorXd, 5> K;
  int d = y.size();
  Eigen::MatrixXd Coeffs = ButcherMatrix();
  for (int i = 0; i < 5; i++) {
    // Calculate coefficients from previous stages:
    Eigen::VectorXd b = Eigen::VectorXd::Zero(d);
    for (int j = 0; j < i; j++) {
      b += Coeffs(i, j) * K.at(j);
    }
    b *= h;
    // Compute increment i and store in K.at(i):
    K.at(i) = solveGenIncrementEquation(f, J, y, b, h, rtol, atol);
  }
  return K;
}
/* SAM_LISTING_END_3 */

/* SAM_LISTING_BEGIN_4 */
template <typename Functor, typename Jacobian>
Eigen::VectorXd discEvolSDIRK(Functor &&f, Jacobian &&J,
                              const Eigen::VectorXd &y, double h,
                              double rtol = 1E-6, double atol = 1E-8) {
  // The b weights are in the last row of Coeffs.
  Eigen::MatrixXd Coeffs = ButcherMatrix();
  Eigen::VectorXd Psi;
  // TO DO (0-2.e)
  //====================
  // Your code goes here
  //====================
  return Psi;
}
/* SAM_LISTING_END_4 */

std::vector<Eigen::VectorXd> solveGradientFlow(const Eigen::VectorXd &d,
                                               double lambda,
                                               const Eigen::VectorXd &y,
                                               double T, unsigned int N);

}  // namespace GradientFlow

#endif  // #ifndef GRADIENTFLOW_H_
