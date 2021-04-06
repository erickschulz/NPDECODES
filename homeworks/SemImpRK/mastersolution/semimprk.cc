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

#include "polyfit.h"

namespace SemImpRK {

/* SAM_LISTING_BEGIN_0 */
double cvgRosenbrock() {
  double cvgRate = 0.0;
  // TO DO: (13-2.d) Use polyfit() to estimate the rate of convergence
  // for solveRosenbrock().
  // Final time
  const double T = 10.0;
  // Mesh sizes h=2^{-k} for k in K.
  const Eigen::ArrayXd K = Eigen::ArrayXd::LinSpaced(7, 4, 10);

  // Initial data
  Eigen::Vector2d y0(1., 1.);
  // Parameter and useful matrix for f
  const double lambda = 1;
  Eigen::Matrix2d R;
  R << 0.0, -1.0, 1.0, 0.0;

  // Function and its Jacobian
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

  // Reference mesh size
  const int N_ref = 10 * std::pow(2, 12);
  // Reference solution
  std::vector<Eigen::VectorXd> solref = solveRosenbrock(f, df, y0, N_ref, T);

  Eigen::ArrayXd Error(K.size());
  std::cout << std::setw(15) << "N" << std::setw(16) << "maxerr\n";
  // Main loop: loop over all meshes
  for (unsigned int i = 0; i < K.size(); ++i) {
    // h = 2^{-k} => N = T*h = T*2^k
    int N = T * std::pow(2, K[i]);
    // Get solution
    std::vector<Eigen::VectorXd> sol = solveRosenbrock(f, df, y0, N, T);
    // Compute error
    double maxerr = 0;
    for (unsigned int j = 0; j < sol.size(); ++j) {
      maxerr =
          std::max(maxerr, (sol.at(j) - solref.at((j * N_ref) / N)).norm());
    }

    Error[i] = maxerr;
    std::cout << std::setw(15) << N << std::setw(16) << maxerr << std::endl;
  }
  // Use log(N)=log(T*2^k)=log(T)+log(2)*k to get natural logarithm of N.
  cvgRate = -polyfit(std::log(2.0) * K, Error.log(), 1)(0);
  return cvgRate;
}
/* SAM_LISTING_END_0 */

}  // namespace SemImpRK
