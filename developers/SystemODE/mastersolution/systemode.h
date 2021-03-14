#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

#include "polyfit.h"
#include "rkintegrator.h"

namespace SystemODE {

/* SAM_LISTING_BEGIN_1 */
template <class Function, class State>
void rk4step(Function &&odefun, double h, const State &y0, State &y1) {
  // TO DO: (12-4.c) Implement a single step of RK4 for
  // the ODE y' = odefun(y), starting from y0 and using
  // the step size h. Save the result in y1.
#if SOLUTION
  auto k1 = odefun(y0);
  auto k2 = odefun(y0 + h / 2 * k1);
  auto k3 = odefun(y0 + h / 2 * k2);
  auto k4 = odefun(y0 + h * k3);

  y1 = y0 + h / 6 * k1 + h / 3 * k2 + h / 3 * k3 + h / 6 * k4;
#else
  //====================
  // Your code goes here
  //====================
#endif
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_0 */
double testcvgRK4() {
  // Parameters and initial values:
  double T = 1;
  int n = 5;
  Eigen::VectorXd y0(2 * n);
  for (int i = 0; i < n; ++i) {
    y0(i) = (i + 1.) / n;
    y0(i + n) = -1;
  }

  // TO DO: (12-4.d) Implement the function f from (12-4.a) as a lambda
  // function, i.e. the right hand side of the system of first order ODEs.
#if SOLUTION
  // Define the relevant tridiagonal matrix.
  Eigen::SparseMatrix<double> C(n, n);
  // Each column has at most 3 non-zero entries.
  C.reserve(Eigen::VectorXi::Constant(n, 3));
  C.insert(0, 0) = 2;
  for (int i = 1; i < n; ++i) {
    C.insert(i, i) = 2;
    C.insert(i, i - 1) = -1;
    C.insert(i - 1, i) = -1;
  }
  C.makeCompressed();

  auto f = [n, C](Eigen::VectorXd y) {
    // The system of ODEs is y' = f(y) with y = [u;v],
    // and f(y) = f([u;v]) = [v;C^{-1}g(u)]
    Eigen::VectorXd fy(2 * n);
    fy.head(n) = y.tail(n);

    Eigen::VectorXd g(n);
    g(0) = y(0) * (y(1) + y(0));
    g(n - 1) = y(n - 1) * (y(n - 1) + y(n - 2));
    for (int i = 1; i < n - 1; ++i) {
      g(i) = y(i) * (y(i - 1) + y(i + 1));
    }

    Eigen::SparseLU<Eigen::SparseMatrix<double>> Csolver;
    Csolver.compute(C);
    fy.tail(n) = Csolver.solve(g);
    return fy;
  };
#else
  //====================
  // Your code goes here
  //====================
#endif

  double conv_rate = 0;
  // Table header
  std::cout << std::setw(8) << "N" << std::setw(20) << "Error" << std::endl;

  // TO DO: (12-4.d) Tabulate the error at time T=1, using
  // N=2,4,...,1024 RK4 steps. Use $N=2^{12}$ steps to calculate
  // the "exact" solution. Then, estimate the algebraic convergence
  // rate of the errors.
  // HINT: You can use polyfit() to calculate the convergence rate.
#if SOLUTION
  // For reference, calculate solution using $N=2^{12}$ steps.
  int N_exact = std::pow(2, 12);
  double h = T / N_exact;  // step size
  Eigen::VectorXd yT_exact = y0;
  for (int step = 0; step < N_exact; step++) {
    // yT_exact is the solution at time t=h*step
    Eigen::VectorXd y_tmp;
    rk4step<std::function<Eigen::VectorXd(Eigen::VectorXd)>, Eigen::VectorXd>(
        f, h, yT_exact, y_tmp);
    yT_exact = y_tmp;
  }

  // Calculate solution using $N=2,...,2^{kmax}$ steps.
  int kmax = 10;
  Eigen::VectorXd Error(kmax);
  for (int k = 0; k < kmax; k++) {
    int N = std::pow(2, k + 1);  // number of steps
    double h = T / N;            // step size
    Eigen::VectorXd yT = y0;
    // Take N RK4 steps:
    for (int step = 0; step < N; step++) {
      // yT is the solution at time t=h*step
      Eigen::VectorXd y_tmp;
      rk4step<std::function<Eigen::VectorXd(Eigen::VectorXd)>, Eigen::VectorXd>(
          f, h, yT, y_tmp);
      yT = y_tmp;
    }
    Error(k) = (yT - yT_exact).norm();
    std::cout << std::setw(8) << N << std::setw(20) << Error(k) << std::endl;
  }

  // Estimate convergence rate
  // Get natural logarithm of N by log(N) = log(2)*log2(N).
  Eigen::VectorXd logN =
      std::log(2) * Eigen::VectorXd::LinSpaced(kmax, 1, kmax);
  Eigen::VectorXd coeffs = polyfit(logN, Error.array().log(), 1);
  conv_rate = -coeffs(0);
#else
  //====================
  // Your code goes here
  //====================
#endif

  return conv_rate;
}
/* SAM_LISTING_END_0 */

}  // namespace SystemODE
