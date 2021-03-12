#include <Eigen/Dense>
#include <Eigen/LU>
#include <iostream>
#include <utility>

#include "initcondlv.h"

int main() {
  /* SAM_LISTING_BEGIN_2 */
  // Initial guess
  Eigen::Vector2d y;
  y << 3, 2;
  // Period we require
  double T = 5;

#if SOLUTION
  // Compute Phi and W from first guess
  std::pair<Eigen::Vector2d, Eigen::Matrix2d> PaW =
      InitCondLV::PhiAndW(y(0), y(1), T);
  // Check error
  Eigen::Vector2d F = PaW.first - y;
  Eigen::Matrix2d DF;

  // Until we are happy with our approximation
  while (F.norm() > 1e-5) {
    // Test current guess
    PaW = InitCondLV::PhiAndW(y(0), y(1), T);
    // Find out error (we want to find a zero of the error)
    F = PaW.first - y;
    // Find out Jacobian
    DF = PaW.second - Eigen::MatrixXd::Identity(2, 2);
    // Use newton iteration
    y = y - DF.lu().solve(F);
  }

  std::cout << "The obtained initial condition is: " << std::endl
            << y << std::endl;
  PaW = InitCondLV::PhiAndW(y(0), y(1), 100);

  std::cout << "y(100) = " << std::endl << PaW.first << std::endl;
#else   // TEMPLATE
  // TODO: Apply the Newton method to find initial data
  // giving solutions with period equal to 5.
#endif  // TEMPLATE
  /* SAM_LISTING_END_2 */

  return 0;
}
