#include <Eigen/Core>
#include <iostream>
#include <vector>

#include "crossprod.h"

int main() {
  double T = 1.;
  int N = 1;

  Eigen::Vector3d y0(0.1, 0.2, 0.4);

  auto f = [](Eigen::Vector3d y) -> Eigen::Vector3d {
    return Eigen::Vector3d(y(0) * y(1), y(1) * y(2), y(2) - y(0));
  };

  auto Jf = [](Eigen::Vector3d y) -> Eigen::Matrix3d {
    Eigen::Matrix3d J;
    J << y(1), y(0), 0, 0, y(2), y(1), -1, 0, 1;
    return J;
  };
  // test implicit midpoint
  std::vector<Eigen::VectorXd> test_imp =
      CrossProd::solve_imp_mid(f, Jf, T, y0, N);
  std::cout << "Implicit midpoint:\n"
            << test_imp.back() << std::endl
            << std::endl;

  // test linear implicit midpoint
  std::vector<Eigen::VectorXd> test_lin =
      CrossProd::solve_lin_mid(f, Jf, T, y0, N);
  std::cout << "Implicit linear midpoint:\n"
            << test_lin.back() << std::endl
            << std::endl;

  // CrossProd::tab_crossprod();

  return 0;
}
