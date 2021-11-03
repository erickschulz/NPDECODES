/**
 * @file taylorode.cc
 * @brief NPDE homework TaylorODE
 * @author ?, Philippe Peter
 * @date 24.03.2021
 * @copyright Developed at ETH Zurich
 */

#include "taylorode.h"

#include <Eigen/Core>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

namespace TaylorODE {

/* SAM_LISTING_BEGIN_1 */
Eigen::Vector2d PredPreyModel::f(const Eigen::Vector2d& y) const {
  //====================
  // Your code goes here
  //====================
  return y;
}

Eigen::Vector2d PredPreyModel::df(const Eigen::Vector2d& y,
                                  const Eigen::Vector2d& z) const {
  //====================
  // Your code goes here
  //====================
  return y;
}

Eigen::Vector2d PredPreyModel::d2f(const Eigen::Vector2d& y,
                                   const Eigen::Vector2d& z) const {
  //====================
  // Your code goes here
  //====================
  return y;
}

std::vector<Eigen::Vector2d> SolvePredPreyTaylor(const PredPreyModel& model,
                                                 double T,
                                                 const Eigen::Vector2d& y0,
                                                 unsigned int M) {
  std::vector<Eigen::Vector2d> res;
  res.reserve(M + 1);

  //====================
  // Your code goes here
  //====================
  return res;
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
double TestCvgTaylorMethod() {
  // initialize parameters for the model:
  double T = 10;               // final time
  Eigen::Vector2d y0(100, 5);  // initial condition
  Eigen::Vector2d yex(0.319465882659820,
                      9.730809352326228);  // reference solution
  double alpha1 = 3.0;
  double alpha2 = 2.0;
  double beta1 = 0.1;
  double beta2 = 0.1;
  PredPreyModel model(alpha1, alpha2, beta1, beta2);

  // Initialize parameters for the convergence study
  unsigned int M0 = 128;    // Minimum number of timesteps
  unsigned int numRef = 8;  // Number of refinements

  // Convergence study
  Eigen::ArrayXd error(numRef);
  Eigen::ArrayXd M(numRef);
  // Run convergence study
  //====================
  // Your code goes here
  //====================

  PrintErrorTable(M, error);

  // Estimate convergence rate based on linear regression
  //====================
  // Your code goes here
  //====================
  return 0.0;
}
/* SAM_LISTING_END_2 */

void PrintErrorTable(const Eigen::ArrayXd& M, const Eigen::ArrayXd& error) {
  std::cout << std::setw(15) << "M" << std::setw(15) << "error" << std::setw(15)
            << "rate" << std::endl;

  for (unsigned int i = 0; i < M.size(); ++i) {
    std::cout << std::setw(15) << M(i) << std::setw(15) << error(i);
    if (i > 0) {
      std::cout << std::setw(15) << std::log2(error(i - 1) / error(i));
    }
    std::cout << std::endl;
  }
}

}  // namespace TaylorODE
