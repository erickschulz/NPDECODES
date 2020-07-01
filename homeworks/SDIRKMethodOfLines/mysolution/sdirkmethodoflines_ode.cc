/** @file
 * @brief NPDE SDIRKMethodOfLines
 * @author Erick Schulz
 * @date 12/04/2019
 * @copyright Developed at ETH Zurich
 */

#include "sdirkmethodoflines_ode.h"

namespace SDIRKMethodOfLines {

/* SAM_LISTING_BEGIN_1 */
std::vector<double> sdirk2SteppingLinScalODE(unsigned int m) {
  std::vector<double> sol_vec;
  double step_size = 2.0 / m;  // "tau" in this example
  // Initial conditions
  sol_vec.push_back(1.0);
// Discrete evolution operator. For the SDIRK-2 method applied to the
// scalar linear ODE (d/dt)y = -y, this turns out to be a scalar valued
// function depending on the step_size. Since we take equidistant step sizes
// in this example, the action of the evolution operator is simply
// multiplication by a constant double
  //=============================================
  // Your task is to modify the follow loop:
  for (int i = 1; i < m + 1; i++) {
    sol_vec.push_back(0.0);
  }
//=============================================
  return sol_vec;
}
/* SAM_LISTING_END_1 */

void sdirk2ScalarODECvTest() {
  int nIter = 10;  // total number of iterations
  unsigned int m;  // number of equidistant steps

  // Evaluating the maximal error (discrete infinity norm) between the approx
  // solutions as given by SDIRK-2 and the exact solution vector computed from
  // the anlytic formula
  /* SAM_LISTING_BEGIN_2 */
  double diff;  // temporary variable used to compute error at various nodes
  double max_norm_errors[nIter];  // errors vector for all approx. sols
  std::vector<double> approx_sol_vec;
  for (int k = 0; k < nIter; k++) {
    m = 10 * std::pow(2, k);

    // Creating exact solution vector. This vector is created by evaluating the
    // exact solution using the analytic formula y(t) = exp(-t) at the
    // equidistant nodes of the time steps.
    double step_size = 2.0 / m;  // "tau" in this example
    std::vector<double> exact_sol_vec;
    for (int i = 0; i < m + 1; i++) {
      exact_sol_vec.push_back(std::exp(-i * step_size));
    }
    // Computing approximate solution
    approx_sol_vec = sdirk2SteppingLinScalODE(m);
    // Computing the error in the maximum norm
    for (int i = 0; i < m + 1; i++) {
      max_norm_errors[k] = 0;
      diff = std::abs(approx_sol_vec.at(i) - exact_sol_vec.at(i));
      max_norm_errors[k] =
          (diff > max_norm_errors[k]) ? diff : max_norm_errors[k];
    }
  }
  // Computing rates of convergence
  double rates[nIter - 1];
  double avg_rate = 0.0;
  for (int k = 0; k < nIter - 1; k++) {
    rates[k] = log2(max_norm_errors[k] / max_norm_errors[k + 1]);
    avg_rate += rates[k];
  }
  avg_rate = avg_rate / (nIter - 1);
  /* SAM_LISTING_END_2 */

  // Printing results
  std::cout << "\n" << std::endl;
  std::cout << "*********************************************************"
            << std::endl;
  std::cout << "         Convergence of SDIRK-2 Method           " << std::endl;
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

}  // namespace SDIRKMethodOfLines
