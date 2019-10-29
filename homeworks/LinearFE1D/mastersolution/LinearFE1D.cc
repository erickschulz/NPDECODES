/**
 * @ file LinearFE1D.cc
 * @ brief NPDE homework LinearFE1D code
 * @ author Christian Mitsch
 * @ date 01.03.2019
 * @ copyright Developed at ETH Zurich
 */

#include "solve_LinearFE1D.h"

int main() {
  // There is no main function to be implemented in this exercise but feel free
  // to use this main to call and test your function

  // BEGIN_SOLUTION
  // Vector mesh = Vector::LinSpaced(11, 0., 1.);
  Vector mesh(9);
  mesh << 0.0, 0.12, 0.2, 0.25, 0.5, 0.7, 0.79, 0.80, 1.0;
  auto alpha = [](double x) { return x; };
  auto f = [](double x) { return x; };
  auto gamma = [](double x) { return x; };
  Eigen::VectorXd uA, uB, uC;
  uA = LinearFE1D::solveA(mesh, gamma, f);
  uB = LinearFE1D::solveB(mesh, alpha, f, 0.1, 0.5);
  uC = LinearFE1D::solveC(mesh, alpha, gamma);

  std::cout << "solveA:" << std::endl;
  std::cout << uA << std::endl;
  std::cout << "solveB:" << std::endl;
  std::cout << uB << std::endl;
  std::cout << "solveC:" << std::endl;
  std::cout << uC << std::endl;

  // END_SOLUTION

  std::cout
      << "You can use this main file to call the function LinearFE1D::solveB"
      << std::endl;
}
