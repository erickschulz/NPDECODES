/**
 * @ file semilagrangian_main.cc
 * @ brief NPDE homework TEMPLATE MAIN FILE
 * @ author
 * @ date May 2022
 * @ copyright Developed at SAM, ETH Zurich
 */

#include <Eigen/Core>
#include <iostream>

#include "semilagrangian.h"

int main() {

  SemiLagrangian::testfloor_and_division();
  //double T = M_PI / 2.;
  double T = 0.5;
  //auto v = [](const Eigen::Vector2d& x) {
  //  return (Eigen::Vector2d() << -(x(1) - 0.5), x(0) - 0.5).finished();
  //};
   auto v = [](const Eigen::Vector2d& x){
    return Eigen::Vector2d(-x(1),x(0));
  };
  
  auto u0 = [](const Eigen::Vector2d& x) {
    Eigen::Vector2d x0 = x;
    x0(0) -= 0.25;
    x0(1) -= 0.5;
    if (x0.norm() < 0.25) {
      return std::pow(std::cos(2. * M_PI * x0.norm()), 2);
    } else {
      return 0.;
    }
  };

  std::cout << "M"
            << "\t"
            << "K"
            << "\t"
            << "error" << std::endl;

  for (int M = 10; M <= 640; M *= 2) {
    int N = (M - 1) * (M - 1);
    Eigen::MatrixXd grid = SemiLagrangian::findGrid(M);

    for (int K = 10; K <= 640; K *= 2) {
      Eigen::VectorXd u = SemiLagrangian::semiLagrangePureTransport(M, K, T);

      Eigen::VectorXd u_ex(N);
      for (int i = 0; i < grid.cols(); ++i) {
        Eigen::Vector2d x = grid.col(i);
        u_ex(i) = SemiLagrangian::solveTransport(x, 640, T, v, u0);
      }

      double err = (u - u_ex).cwiseAbs().maxCoeff();
      std::cout << M << "\t" << K << "\t" << err << std::endl;
    }
  }
}