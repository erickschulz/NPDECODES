#include <Eigen/Core>
#include <iostream>

#include "nlmatode.h"

int main() {
  double T = 1;
  unsigned int n = 3;

  Eigen::MatrixXd Y0(n, n);
  Y0 << 1, 1, 0, 0, 3, 2, 1, 5, 2;

  Eigen::MatrixXd YT = NLMatODE::matode(Y0, T);

  std::cout << YT << std::endl << std::endl;
  Eigen::MatrixXd M0(n, n);
  M0 << 1, 1, 0, 1, 3, 1, 0, 1, 1;
  // check whether invariant is perserved or not
  bool is_invariant = NLMatODE::checkinvariant(Y0, T);
  bool is_invar_sym = NLMatODE::checkinvariant(M0, T);
  std::cout << "Test whether invariant was preserved or not..." << std::endl;

  if (is_invariant) {
    std::cout << "Invariant for Y0 preserved." << std::endl;
  } else {
    std::cout << "Invariant for Y0 NOT preserved." << std::endl;
  }
  if (is_invar_sym) {
    std::cout << "Invariant for M0 preserved." << std::endl;
  } else {
    std::cout << "Invariant for M0 NOT preserved." << std::endl;
  }

  double rate = NLMatODE::cvgDiscreteGradientMethod();
  std::cout << "\nThe fitted rate for the discrete gradient method is:\n"
            << rate << std::endl;
  return 0;
}
