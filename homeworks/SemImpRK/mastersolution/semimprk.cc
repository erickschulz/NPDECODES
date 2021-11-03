/**
 * @file semimprk.cc
 * @brief NPDE homework SemImpRK code
 * @author Unknown, Oliver Rietmann
 * @date 04.04.2021
 * @copyright Developed at ETH Zurich
 */

#include "semimprk.h"

#include <Eigen/Core>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

#include "../../../lecturecodes/helperfiles/polyfit.h"

namespace SemImpRK {

/* SAM_LISTING_BEGIN_0 */
double CvgRosenbrock() {
  double cvgRate = 0.0;
  // Use polyfit to estimate the rate of convergence
  // for SolveRosenbrock.
  // Final time
  const double T = 10.0;
  // Timestep sizes h=2^{-k} for which discretization erorrs should be computed
  const Eigen::ArrayXd K = Eigen::ArrayXd::LinSpaced(7, 4, 10);
  // Initial state vector
  Eigen::Vector2d y0(1., 1.);

  // Parameter and useful matrix for f
  const double lambda = 1;
  Eigen::Matrix2d R;
  R << 0.0, -1.0, 1.0, 0.0;
  // Lambda functions giving the right-hand side function f
  // and its Jacobian df
  auto f = [&R, &lambda](Eigen::Vector2d y) {
    return R * y + lambda * (1.0 - y.squaredNorm()) * y;
  };
  auto df = [&lambda](Eigen::Vector2d y) {
    double x = 1 - y.squaredNorm();
    Eigen::Matrix2d J;
    J << lambda * x - 2 * lambda * y(0) * y(0), -1 - 2 * lambda * y(1) * y(0),
        1 - 2 * lambda * y(1) * y(0), lambda * x - 2 * lambda * y(1) * y(1);
    return J;
  };

  // timestep size h=2^{-12} (=> M = T*2^{12}) used to compute a reference
  // solution
  const int M_ref = T * std::pow(2, 12);
  // Reference solution
  std::vector<Eigen::VectorXd> solref = SolveRosenbrock(f, df, y0, M_ref, T);
  // Arrays for recording error norms and numbers of timesteps
  Eigen::ArrayXd Error(K.size());
  Eigen::ArrayXd M(K.size());
  // Loop over all temporal meshes
  for (unsigned int i = 0; i < K.size(); ++i) {
    // h = 2^{-k} => M = T*h = T*2^k
    M(i) = T * std::pow(2, K[i]);

    // Compute solution
    std::vector<Eigen::VectorXd> sol = SolveRosenbrock(f, df, y0, M(i), T);

    // Compute error
    double maxerr = 0;
    for (unsigned int j = 0; j < sol.size(); ++j) {
      maxerr = std::max(maxerr, (sol[j] - solref[j * M_ref / M(i)]).norm());
    }
    Error(i) = maxerr;
  }
  // Print Error table
  std::cout << std::setw(15) << "M" << std::setw(16) << "maxerr \n";
  for (int i = 0; i < M.size(); ++i) {
    std::cout << std::setw(15) << M(i) << std::setw(15) << Error(i)
              << std::endl;
  }
  // Estimate convergence rate using linear regression provided by the auxiliary
  // function polyfit()
  cvgRate = -polyfit(M.log(), Error.log(), 1)(0);
  return cvgRate;
}
/* SAM_LISTING_END_0 */

}  // namespace SemImpRK
