/**
 * @file radauthreetimesteppingode.cc
 * @brief NPDE homework RadauThreeTimestepping
 * @author Erick Schulz
 * @date 08/04/2019
 * @copyright Developed at ETH Zurich
 */

#include "radauthreetimesteppingode.h"

#include <cmath>
#include <iostream>
#include <vector>

namespace RadauThreeTimestepping {

/* SAM_LISTING_BEGIN_1 */
std::vector<double> twoStageRadauTimesteppingLinScalODE(unsigned int m) {
  std::vector<double> sol_vec;
#if SOLUTION
  double step_size = 5.0 / m;  // Timestep "tau"
  sol_vec.push_back(1.0);      // Initial value

  // Discrete evolution operator: For the two-stage Radau method applied to the
  // scalar linear ODE (d/dt)y = -y, this turns out to be a scalar valued
  // function depending on the step_size. Since we take equidistant step sizes
  // in this example, the action of the evolution operator is simply
  // multiplication by a constant double
  double evolution_op =
      (1 - (step_size * (1 + step_size / 6.0)) /
               ((1 + step_size * 5.0 / 12.0) * (1 + step_size / 4.0) +
                step_size * step_size / 16.0));
  // Compute discrete evolution by applying the evolution operator at each step
  for (int i = 1; i < m + 1; i++) {
    sol_vec.push_back(evolution_op * sol_vec.at(i - 1));
  }
#else
  //====================
  // Your code goes here
  //====================
#endif
  return sol_vec;
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
void testConvergenceTwoStageRadauLinScalODE() {
  constexpr int nIter = 10;       // total number of iterations
  double max_norm_errors[nIter];  // errors vector for all approx. sols
  double rates[nIter - 1];        // The rates of convergence
  double avg_rate = 0.0;  // The average rate of convergence over all iterations

#if SOLUTION
  // Error between the approx solutions as given by the two stage Radau method
  // and the exact solution vector computed from the anlytic formula vector
  // computed from the anlytic formula
  double diff;  // temporary variable used to compute error at various nodes
  std::vector<double> approx_sol_vec;
  for (int k = 0; k < nIter; k++) {
    unsigned int m = 10 * std::pow(2, k);  // number of equidistant steps
    double step_size = 5.0 / m;            // time step `tau`
    // Creating exact solution vector. This vector is created by evaluating the
    // exact solution using the analytic formula y(t) = exp(-t) at the
    // equidistant nodes of the time steps.
    std::vector<double> exact_sol_vec;
    for (int i = 0; i < m + 1; i++) {
      exact_sol_vec.push_back(std::exp(-i * step_size));
    }
    // Computing approximate solution
    approx_sol_vec = twoStageRadauTimesteppingLinScalODE(m);
    // Computing the error in the maximum norm
    for (int i = 0; i < m + 1; i++) {
      max_norm_errors[k] = 0;
      diff = std::abs(approx_sol_vec.at(i) - exact_sol_vec.at(i));
      max_norm_errors[k] =
          (diff > max_norm_errors[k]) ? diff : max_norm_errors[k];
    }
  }
  // Computing rates of convergence
  for (int k = 0; k < nIter - 1; k++) {
    rates[k] = log2(max_norm_errors[k] / max_norm_errors[k + 1]);
    avg_rate += rates[k];
  }
  avg_rate = avg_rate / (nIter - 1);
#else
  //====================
  // Your code goes here
  //====================
#endif
  /* SAM_LISTING_END_2 */

  // Printing results
  std::cout << "\n" << std::endl;
  std::cout << "*********************************************************"
            << std::endl;
  std::cout << "         Convergence of two-stage Radau Method           "
            << std::endl;
  std::cout << "*********************************************************"
            << std::endl;
  std::cout << "--------------------- RESULTS ---------------------------"
            << std::endl;
  std::cout << "Iteration"
            << "\t| Nsteps"
            << "\t| error"
            << "\t\t| rates" << std::endl;
  std::cout << "---------------------------------------------------------"
            << std::endl;
  for (int k = 0; k < nIter; k++) {
    std::cout << k << "\t"
              << "\t|" << 10 * std::pow(2, k) << "\t\t|" << max_norm_errors[k];
    if (k > 0) {
      std::cout << "\t|" << rates[k - 1];
    }
    std::cout << "\n";
  }
  std::cout << "---------------------------------------------------------"
            << std::endl;
  std::cout << "Average rate of convergence: " << avg_rate << "\n" << std::endl;
}

}  // namespace RadauThreeTimestepping
