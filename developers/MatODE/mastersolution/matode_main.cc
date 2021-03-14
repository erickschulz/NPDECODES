#include <Eigen/Core>
#include <iostream>

#include "matode.h"

int main() {
  MatODE::checkOrthogonality();

  double T = 1.0;
  unsigned int n = 3;
  Eigen::MatrixXd M(n, n);
  M << 8.0, 1.0, 6.0, 3.0, 5.0, 7.0, 4.0, 9.0, 2.0;

  std::cout << "Test implementation of ode45" << std::endl;

  std::cout << "M = " << std::endl << M << std::endl;
  Eigen::MatrixXd N = MatODE::matode(M, T);
  std::cout << "N = " << std::endl << N << std::endl;

  std::cout << "Test whether invariant was preserved or not" << std::endl;
  bool is_invariant = MatODE::checkinvariant(N, T);

  if (is_invariant) {
    std::cout << "Invariant was preserved." << std::endl;
  } else {
    std::cout << "Invariant was NOT preserved." << std::endl;
  }
}
