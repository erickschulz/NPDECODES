#include <Eigen/Core>
#include <cmath>
#include <functional>
#include <iostream>

#include "systemode.h"

int main() {
  Eigen::VectorXd y0(3), y1(3);
  y0 << -1, 1, 2;
  double h = 0.5;
  auto f = [](Eigen::VectorXd y) {
    Eigen::VectorXd fy(3);
    fy << y(1) * y(2), y(0) * y(1), 3 * y(2);
    return fy;
  };
  SystemODE::rk4step<std::function<Eigen::VectorXd(Eigen::VectorXd)>,
                     Eigen::VectorXd>(f, h, y0, y1);
  std::cout << "Test of rk4step() gives y1 = " << y1.transpose() << std::endl;

  std::cout << "Convergence rate: "
            << std::round(std::abs(SystemODE::testcvgRK4())) << std::endl;

  return 0;
}
