#ifndef STABRK3_H_
#define STABRK3_H_

/**
 * @file stabrk3.h
 * @brief NPDE homework StabRK3 code
 * @author Oliver Rietmann, Philippe Peter
 * @date 13.04.2021
 * @copyright Developed at ETH Zurich
 */

#include <Eigen/Core>
#include <vector>

namespace StabRK3 {

/**
 * @brief Solve the predator-prey IVP using the RK-SSM and uniform timesteps.
 * @param y0 initial condition
 * @param T final time
 * @param N number of uniform timesteps
 * @return approximate solution at final time T
 */
Eigen::Vector2d PredPrey(Eigen::Vector2d y0, double T, unsigned N);

/**
 * @brief Study the asymptotic convergence behavior of the method.
 */
void SimulatePredPrey();

/**
 * @brief Helper function, prints an error table.
 */
void PrintErrorTable(const Eigen::ArrayXd& N, const Eigen::ArrayXd& error);

}  // namespace StabRK3

#endif  // #ifndef STABRK3_H_
