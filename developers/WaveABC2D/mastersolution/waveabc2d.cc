/** @file
 * @brief NPDE WaveABC2D
 * @author Erick Schulz
 * @date 11/12/2019
 * @copyright Developed at ETH Zurich
 */

#include "waveabc2d.h"

#include <cmath>
#include <iomanip>
#include <iostream>
// Eigen includes
#include <Eigen/Core>

namespace WaveABC2D {

/** @brief The evolution is performed by successively solving a linear sytem
 *
 *                        A * x(k) = b(k-1)
 *
 * at each time step. The matrix A is independent of time and completely
 * determined by the equation's parameter epsilon and the fixed uniform step
 * size of the method. Here, x will be a (M + 1) x 2 matrix storing the discrete
 * solution in its first column and the discrete time derivative in the second
 * column. The discrete time derivative is needed, because it is introduced as
 * an unknown in the implicit method to reduce the order of the ODE.
 *
 * @param epsilon Equation parameter in [0,1].
 * @param M Number of timesteps.*/
/* SAM_LISTING_BEGIN_1 */
Eigen::VectorXd scalarImplicitTimestepping(double epsilon, unsigned int M) {
  /* PROBLEM SETUP */
  double step_size = 1.0 / M;  // time step "tau"
  // Rows contain sequence of solutions states
  Eigen::MatrixXd x(M + 1, 2);
  // INITIAL DATA
  x(0, 0) = 1.0;
  x(0, 1) = 1.0;

#if SOLUTION
  Eigen::MatrixXd A(2, 2);
  Eigen::VectorXd b(2);
  // clang-format off
    A << epsilon + 0.5*step_size*(1-epsilon), 0.5*step_size,
         -0.5*step_size,                      1.0;
  // clang-format on
  // TIMESTEPPING
  for (int k = 1; k < M + 1; k++) {
    // Update right hand side
    b(0) = (epsilon - 0.5 * step_size * (1 - epsilon)) * x(k - 1, 0) -
           0.5 * step_size * x(k - 1, 1);
    b(1) = x(k - 1, 1) + 0.5 * step_size * x(k - 1, 0);
    // Solve the linear system
    x.row(k) = A.fullPivLu().solve(b);
  }
#else
//====================
// Your code goes here
//====================
#endif

  return x.col(1);
}  // scalarImplicitTimestepping
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
void testConvergenceScalarImplicitTimestepping() {
  std::cout << "Testing convergence of the implicit method." << std::endl;
  int nIter = 13;    // total number of iterations
  unsigned int M;    // number of equidistant steps
  double step_size;  // time step "tau"
  double epsilon = 0.2;

  // Error between the approx solutions as given by the implicit method
  // and the exact solution vector computed from the anlytic formula vector
  // computed from the anlytic formula
  std::cout << "Computing approximate solutions..." << std::endl;
  double diff;  // temporary variable used to compute error at various nodes
  double max_norm_errors[nIter];  // errors vector for all approx. sols
  Eigen::VectorXd approx_sol_vec;
  Eigen::VectorXd exact_sol_vec(M + 1);
  unsigned int M_stored[13] = {10,  20,  30,  40,  50,  60, 80,
                               100, 160, 200, 320, 500, 640};
  for (int k = 0; k < nIter; k++) {
    // M = 10 * std::pow(2, k);
    M = M_stored[k];
    step_size = 1.0 / M;
    // Creating exact solution vector for epsilon = 1/5. This vector is created
    // by evaluating the exact solution using the analytic formula
    //                      x(t) = exp(-2*t)*cos(t)
    // at the equidistant nodes of the time steps.
    for (int i = 0; i < M + 1; i++) {
      exact_sol_vec(i) =
          std::exp(-2 * i * step_size) *
          (3 * std::sin(i * step_size) + std::cos(i * step_size));
    }
    // Computing approximate solution
    approx_sol_vec = scalarImplicitTimestepping(epsilon, M);
    // Computing the error in the maximum norm
    for (int i = 0; i < M + 1; i++) {
      max_norm_errors[k] = 0;
      diff = std::abs(approx_sol_vec(i) - exact_sol_vec(i));
      max_norm_errors[k] =
          (diff > max_norm_errors[k]) ? diff : max_norm_errors[k];
    }
  }
  // Computing rates of convergence
  /*std::cout << "Computing convergence rates..." << std::endl;
    double rates[nIter - 1];
    double avg_rate = 0.0;
    for (int k = 0; k < nIter - 1; k++) {
      rates[k] = log2(max_norm_errors[k] / max_norm_errors[k + 1]);
      avg_rate += rates[k];
    }
    avg_rate = avg_rate / (nIter - 1);
  */
  /* SAM_LISTING_END_2 */

  // Printing results
  std::cout << "\n" << std::endl;
  std::cout << "*********************************************************"
            << std::endl;
  std::cout << "         Convergence of the Implicit Method              "
            << "\n                  (epsilon = " << epsilon << ")" << std::endl;
  std::cout << "*********************************************************"
            << std::endl;
  std::cout << "--------------------- RESULTS ---------------------------"
            << std::endl;
  std::cout << "Iteration"
            << "\t| Nsteps"
            << "\t| error" << std::endl;
  //            << "\t\t| rates" << std::endl;
  std::cout << "---------------------------------------------------------"
            << std::endl;
  for (int k = 0; k < nIter; k++) {
    std::cout << k << "\t"
              << "\t|" << M_stored[k] /*10 * std::pow(2, k)*/ << "\t\t|"
              << max_norm_errors[k];
    /*  if (k > 0) {
          std::cout << "\t|" << rates[k - 1];
        }
    */
    std::cout << "\n";
  }
  std::cout << "---------------------------------------------------------"
            << std::endl;
  /*    std::cout << "Average rate of convergence: " << avg_rate << "\n" <<
   * std::endl;
   */

}  // testConvergenceScalarImplicitTimestepping
/* SAM_LISTING_END_2 */

/* Implementing member functions of class progress_bar */
void progress_bar::write(double fraction) {
  // clamp fraction to valid range [0,1]
  if (fraction < 0)
    fraction = 0;
  else if (fraction > 1)
    fraction = 1;

  auto width = bar_width - message.size();
  auto offset = bar_width - static_cast<unsigned>(width * fraction);

  os << '\r' << message;
  os.write(full_bar.data() + offset, width);
  os << " [" << std::setw(3) << static_cast<int>(100 * fraction) << "%] "
     << std::flush;
}

}  // namespace WaveABC2D
