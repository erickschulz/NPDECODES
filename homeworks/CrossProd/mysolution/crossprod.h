#ifndef CROSSPROD_H_
#define CROSSPROD_H_

/**
 * @file crossprod.h
 * @brief NPDE homework CrossProd code
 * @author Unknown, Oliver Rietmann
 * @date 31.03.2021
 * @copyright Developed at ETH Zurich
 */

#include <Eigen/Core>
#include <Eigen/LU>
#include <utility>
#include <vector>

#include "implicitrkintegrator.h"

namespace CrossProd {

//! \tparam Function type for function implementing the rhs function.
//! Must have Eigen::VectorXd operator()(Eigen::VectorXd x)
//! \tparam Jacobian type for function implementing the Jacobian of f.
//! Must have Eigen::MatrixXd operator()(Eigen::VectorXd x)
//! \param[in] f function handle for rhs in y' = f(y), e.g. implemented
//! using lambda function.
//! \param[in] Jf function handle for Jf, e.g. implemented using lambda function
//! \param[in] T final time T
//! \param[in] y0 initial data y(0) = y0 for y' = f(y)
//! \param[in] N number of steps to perform.
/* SAM_LISTING_BEGIN_1 */
template <class Function, class Jacobian>
std::vector<Eigen::VectorXd> solve_imp_mid(Function &&f, Jacobian &&Jf,
                                           double T, const Eigen::VectorXd &y0,
                                           unsigned int N) {
  std::vector<Eigen::VectorXd> res(N + 1);
  // TO DO (13-1.e): Construct the implicit mid-point method with the class
  // implicit_RKIntegrator and execute the .solve() method.
  // Return the vector containing all steps including initial and final value.
  //====================
  // Your code goes here
  //====================
  return res;
}
/* SAM_LISTING_END_1 */

//! \tparam Function type for function implementing the rhs function.
//! Must have Eigen::VectorXd operator()(Eigen::VectorXd x)
//! \tparam Jacobian type for function implementing the Jacobian of f.
//! Must have Eigen::MatrixXd operator()(Eigen::VectorXd x)
//! \param[in] f function handle for rhs in y' = f(y), e.g. implemented
//! using lambda function.
//! \param[in] Jf function handle for Jf, e.g. implemented using lambda function
//! \param[in] T final time T
//! \param[in] y0 initial data y(0) = y0 for y' = f(y)
//! \param[in] N number of steps to perform.
/* SAM_LISTING_BEGIN_2 */
template <class Function, class Jacobian>
std::vector<Eigen::VectorXd> solve_lin_mid(Function &&f, Jacobian &&Jf,
                                           double T, const Eigen::VectorXd &y0,
                                           unsigned int N) {
  std::vector<Eigen::VectorXd> res(N + 1);
  // TO DO (13-1.g): Implement the linear implicit mid-point method for
  // an autonomous ODE y' = f(y), y(0) = y0. Return the vector containing
  // all steps including initial and final value.
  //====================
  // Your code goes here
  //====================
  return res;
}
/* SAM_LISTING_END_2 */

void tab_crossprod();

}  // namespace CrossProd

#endif  // #define CROSSPROD
