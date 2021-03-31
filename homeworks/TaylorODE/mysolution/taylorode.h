#ifndef TAYLORODE_H
#define TAYLORODE_H
/**
 * @file taylorode.h
 * @brief NPDE homework TaylorODE
 * @author ?, Philippe Peter
 * @date 24.03.2021
 * @copyright Developed at ETH Zurich
 */

#include <Eigen/Dense>
#include <vector>

namespace TaylorODE {

/**
 * @brief Helper class for the evaluation of f, df and d2f of a predator-prey
 * model
 */
class PredPreyModel {
 public:
  /** @brief Construct a predator-prey model based on model parameters */
  PredPreyModel(double alpha1, double alpha2, double beta1, double beta2)
      : alpha1_{alpha1}, alpha2_{alpha2}, beta1_{beta1}, beta2_{beta2} {}

  /** @brief Evaluate f(y) */
  Eigen::Vector2d f(const Eigen::Vector2d& y) const;

  /** @brief Evaluate df(y)*z */
  Eigen::Vector2d df(const Eigen::Vector2d& y, const Eigen::Vector2d& z) const;

  /** @brief Evaluate d2f(y)(z,z) (second derivative bilinear mapping)*/
  Eigen::Vector2d d2f(const Eigen::Vector2d& y, const Eigen::Vector2d& z) const;

 private:
  const double alpha1_;
  const double alpha2_;
  const double beta1_;
  const double beta2_;
};

/**
 * @brief Compute the solution of the predator-prey ode
 * @param model Predator-prey model, allows evaluation of f, df, d2f
 * @param T Final time
 * @param y0 Initial state
 * @param M Number of equidistant timesteps
 * @return Vector containing all values y^m (including initial and final value)
 */
std::vector<Eigen::Vector2d> SolvePredPreyTaylor(const PredPreyModel& model,
                                                 double T,
                                                 const Eigen::Vector2d& y0,
                                                 unsigned int M);

/**
 * @brief Determine the order of convergence for the Taylor expansion method
 * applied to the predator-prey model
 * @return  Empiric rate of convergence determined by linear regression
 */
double TestCvgTaylorMethod();

/**
 * @brief Helper function that prints an error table
 */
void PrintErrorTable(const Eigen::ArrayXd& M, const Eigen::ArrayXd& error);

}  // namespace TaylorODE

#endif  // TAYLORODE_H
