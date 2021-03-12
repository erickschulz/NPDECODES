#include <Eigen/Dense>
#include <cmath>
#include <iostream>

#include "ordnotall.h"

int main() {
  /* SAM_LISTING_BEGIN_2 */
  // Construct data for Butcher schemes
  Eigen::MatrixXd A1 = Eigen::MatrixXd::Zero(1, 1);
  Eigen::VectorXd b1(1);
  b1 << 1;

  Eigen::MatrixXd A2 = Eigen::MatrixXd::Zero(2, 2);
  A2(1, 0) = 1;
  Eigen::VectorXd b2(2);
  b2 << .5, .5;

  Eigen::MatrixXd A3 = Eigen::MatrixXd::Zero(3, 3);
  A3(1, 0) = .5;
  A3(2, 0) = -1;
  A3(2, 1) = 2;
  Eigen::VectorXd b3(3);
  b3 << 1. / 6, 2. / 3, 1. / 6;

  Eigen::MatrixXd A4 = Eigen::MatrixXd::Zero(4, 4);
  A4(1, 0) = .5;
  A4(2, 1) = .5;
  A4(3, 2) = 1;
  Eigen::VectorXd b4(4);
  b4 << 1. / 6, 1. / 3, 1. / 3, 1. / 6;

#if SOLUTION
  // First ODE
  std::cout << std::endl
            << "1. ODE y' = (1-y)y, y(0)=.5" << std::endl
            << std::endl;
  double T = 0.1;
  auto f = [](Eigen::VectorXd y) {
    Eigen::VectorXd fy(1);
    fy << (1. - y(0)) * y(0);
    return fy;
  };
  Eigen::VectorXd y0(1);
  y0 << .5;

  std::cout << "Explicit Euler" << std::endl << std::endl;
  OrdNotAll::errors(f, T, y0, A1, b1);
  std::cout << "Trapezoidal rule" << std::endl << std::endl;
  OrdNotAll::errors(f, T, y0, A2, b2);
  std::cout << "RK order 3" << std::endl << std::endl;
  OrdNotAll::errors(f, T, y0, A3, b3);
  std::cout << "Classical RK order 4" << std::endl << std::endl;
  OrdNotAll::errors(f, T, y0, A4, b4);

  // Second ODE
  std::cout << std::endl
            << "2. ODE y' = |1.1 - y| + 1, y(0)=1" << std::endl
            << std::endl;
  auto f2 = [](Eigen::VectorXd y) {
    Eigen::VectorXd fy(1);
    fy << std::abs(1.1 - y(0)) + 1.;
    return fy;
  };
  y0 << 1;

  std::cout << "Explicit Euler" << std::endl << std::endl;
  OrdNotAll::errors(f2, T, y0, A1, b1);
  std::cout << "Trapezoidal rule" << std::endl << std::endl;
  OrdNotAll::errors(f2, T, y0, A2, b2);
  std::cout << "RK order 3" << std::endl << std::endl;
  OrdNotAll::errors(f2, T, y0, A3, b3);
  std::cout << "Classical RK order 4" << std::endl << std::endl;
  OrdNotAll::errors(f2, T, y0, A4, b4);
#else  // TEMPLATE
  // TODO: call order for all combinations of ODE ans RK method
#endif
  /* SAM_LISTING_END_2 */

  return 0;
}
