/**
 * @file burgersequation.cc
 * @brief NPDE homework BurgersEquation code
 * @author Oliver Rietmann
 * @date 15.04.2019
 * @copyright Developed at ETH Zurich
 */

#include "burgersequation.h"

#include <cmath>

#include <Eigen/Core>

namespace BurgersEquation {

constexpr double PI = 3.14159265358979323846;

double Square(double x) { return x * x; }

double w0(double x) {
  return 0.0 <= x && x <= 1.0 ? Square(std::sin(PI * x)) : 0.0;
}

double f(double x) { return 2.0 / 3.0 * std::sqrt(x * x * x); }

Eigen::VectorXd solveBurgersGodunov(double T, unsigned int N) {
  double h = 5.0 / N;
  double tau = h;
  int m = std::round(T / tau);

  Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(N + 1, -1.0, 4.0);
  Eigen::VectorXd mu = x.unaryExpr(&w0);

  /* SOLUTION_BEGIN */
  /* SOLUTION_END */

  return mu;
}

/**
 * @brief Converts a large vector on  a grid to a smaller vector correponding to
 * a sub-grid.
 *
 * @param mu vector of function values on a spacial grid of size N_large
 * @param N divides the size N_large of mu
 * @return a vector mu_sub of size N, that represents mu on a sub-grid of size N
 */
Eigen::VectorXd reduce(const Eigen::VectorXd &mu, unsigned int N) {
  Eigen::VectorXd mu_sub(N + 1);
  int fraction = mu.size() / N;
  for (int j = 0; j < N + 1; ++j) {
    mu_sub(j) = mu(j * fraction);
  }
  return mu_sub;
}

Eigen::Matrix<double, 3, 4> numexpBurgersGodunov() {
  const unsigned int N_large = 3200;
  Eigen::Vector2d T{0.3, 3.0};
  Eigen::Vector4i N{5 * 10, 5 * 20, 5 * 40, 5 * 80};
  Eigen::Vector4d h;
  for (int i = 0; i < 4; ++i) h(i) = 5.0 / N(i);

  Eigen::Matrix<double, 3, 4> result;
  result.row(0) = h.transpose();

  /* SOLUTION_BEGIN */
  /* SOLUTION_END */

  return result;
}

}  // namespace BurgersEquation
