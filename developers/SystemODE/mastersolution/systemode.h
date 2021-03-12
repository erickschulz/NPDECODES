#include <Eigen/Dense>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

#include "rkintegrator.h"

namespace SystemODE {

/* SAM_LISTING_BEGIN_1 */
template <class Function, class State>
void rk4step(const Function &odefun, double h, const State &y0, State &y1) {
#if SOLUTION
  auto k1 = odefun(y0);
  auto k2 = odefun(y0 + h / 2 * k1);
  auto k3 = odefun(y0 + h / 2 * k2);
  auto k4 = odefun(y0 + h * k3);

  y1 = y0 + h / 6 * k1 + h / 3 * k2 + h / 3 * k3 + h / 6 * k4;
#else   // TEMPLATE
  // TODO: implement a single step of the classical Runge-Kutta method of order
  // 4
#endif  // TEMPLATE
}
/* SAM_LISTING_END_1 */

// This function approximates the order of convergence of the RK scheme defined
// by A and b when applied to the first order system y'=f(y), y(0)=y0. We are
// interested in the error of the solutions at time T.

template <class Function>
void errors(const Function &f, const double &T, const Eigen::VectorXd &y0,
            const MatrixXd &A, const Eigen::VectorXd &b) {
  SystemODE::RKIntegrator<Eigen::VectorXd> rk(A, b);
  std::vector<double> error(15);
  std::vector<double> order(14);
  double sum = 0;
  int count = 0;
  bool test = 1;
  std::vector<Eigen::VectorXd> y_exact = rk.solve(f, T, y0, std::pow(2, 15));

  for (int k = 0; k < 15; k++) {
    int N = pow(2, k + 1);
    std::vector<Eigen::VectorXd> y1 = rk.solve(f, T, y0, N);

    error[k] = (y1[N] - y_exact[pow(2, 15)]).norm();
    std::cout << left << std::setw(3) << std::setfill(' ') << "N = ";
    std::cout << left << std::setw(7) << std::setfill(' ') << N;
    std::cout << left << std::setw(8) << std::setfill(' ') << "Error = ";
    std::cout << left << std::setw(13) << std::setfill(' ') << error[k];

    if (error[k] < y0.size() * 5e-14) test = 0;
    if (k > 0 && test) {
      order[k - 1] = std::log(error[k - 1] / error[k]) / std::log(2.0);
      std::cout << left << std::setw(10) << std::setfill(' ')
                << "Approximated order = " << order[k - 1] << std::endl;
      sum += order[k - 1];
      count = k;
    } else
      std::cout << std::endl;
  }
  std::cout << "Average approximated order = " << sum / count << std::endl
            << std::endl;
}

}  // namespace SystemODE
