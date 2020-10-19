/**
 * @file burgersequation.cc
 * @brief NPDE homework BurgersEquation code
 * @author Oliver Rietmann
 * @date 15.04.2019
 * @copyright Developed at ETH Zurich
 */

#include "burgersequation.h"

#include <Eigen/Core>
#include <cmath>

namespace BurgersEquation {
/* SAM_LISTING_BEGIN_1 */
constexpr double PI = 3.14159265358979323846;

double Square(double x) { return x * x; }

double f(double x) { return 2.0 / 3.0 * std::sqrt(x * x * x); }

Eigen::VectorXd solveBurgersGodunov(double T, unsigned int N) {
  double h = 5.0 / N;           // meshwidth
  double tau = h;               // timestep = meshwidth by CFL condition
  int m = std::round(T / tau);  // no. of timesteps

  // initialize vector with initial nodal values
  Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(N + 1, -1.0, 4.0);
  Eigen::VectorXd mu =
      x.unaryExpr([](double x) {
         return 0.0 <= x && x <= 1.0 ? Square(std::sin(PI * x)) : 0.0;
       }).eval();

  for (int i = 0; i < m; ++i) {
    for (int j = N; 0 < j; --j) {
      // Standard fully discrete evolution based on explicit Euler timestepping
      mu(j) = mu(j) - tau / h * (f(mu(j)) - f(mu(j - 1)));
    }
    // truncation to a finite vector. Only required on one side, because all
    // information flows from left to right.
    mu(0) = 0.0;  // Value of u0 to the left of x=0
  }

  return mu;
}
/* SAM_LISTING_END_1 */

/**
 * @brief Converts a large vector on  a grid to a smaller vector correponding to
 * a sub-grid.
 *
 * @param mu vector of function values on a spacial grid of size N_large
 * @param N divides the size N_large of mu
 * @return a vector mu_sub of size N, that represents mu on a sub-grid of size N
 */
/* SAM_LISTING_BEGIN_2 */
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

  for (int k = 0; k < 2; ++k) {
    Eigen::VectorXd mu_ref = solveBurgersGodunov(T(k), N_large);
    Eigen::Vector4d error;
    for (int i = 0; i < 4; ++i) {
      Eigen::VectorXd mu = solveBurgersGodunov(T(k), N(i));
      Eigen::VectorXd mu_ref_sub = reduce(mu_ref, N(i));
      error(i) = h(i) * (mu - mu_ref_sub).lpNorm<1>();
    }
    result.row(k + 1) = error.transpose();
  }

  return result;
}
/* SAM_LISTING_END_2 */

}  // namespace BurgersEquation
