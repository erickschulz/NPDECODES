#include <Eigen/Core>
#include <iostream>
#include <vector>

#include "odesolve.h"

int main() {
  auto Psi = [](double h, const Eigen::VectorXd& y0) -> Eigen::VectorXd {
    return (1 + h) * y0;
  };

  Eigen::VectorXd y0(1);
  y0 << 1.0;
  double T = 1.0;
  // Test psitilde
  Eigen::VectorXd y1 = ODESolve::psitilde(Psi, 1, 0.1, y0);
  std::cout << "Test psitilde\n" << y1 << std::endl;

  std::cout << "Enter 1 to test equidistant integration\n"
            << "Enter 2 to test adaptive integration\n"
            << "Enter 0 to exit\n";
  unsigned int flag;
  std::cin >> flag;
  switch (flag) {
    case 1: {
      unsigned int N = 8;

      std::cout << "Test equidistant integration" << std::endl;
      std::vector<Eigen::VectorXd> Y1 = ODESolve::odeintequi(Psi, T, y0, N);
      for (int i = 0; i < N; ++i) {
        std::cout << Y1[i](0) << std::endl;
      }
      double rate = ODESolve::testcvpExtrapolatedEuler();
      std::cout << "\nRate = " << rate << std::endl;
      break;
    }

    case 2: {
      std::cout << "\n\nTest adaptive integration" << std::endl;
      std::vector<Eigen::VectorXd> mysol =
          ODESolve::odeintssctrl(Psi, T, y0, 0.01, 1, 10e-5, 10e-5, 10e-5);
      for (int i = 0; i < 8; ++i) {
        std::cout << mysol[i](0) << std::endl;
      }

      ODESolve::solveTangentIVP();
      break;
    }
    default:
      break;
  }

  return 0;
}
