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

// Solve the predator-prey IVP using the RK-SSM and uniform timesteps.
Eigen::Vector2d PredPrey(Eigen::Vector2d y0, double T, unsigned N);

// Study the asymptotic convergence behavior of the method.
void SimulatePredPrey();

// Helper function, prints an error table.
void PrintErrorTable(const Eigen::ArrayXd& N, const Eigen::ArrayXd& error);

}  // namespace StabRK3

#endif  // #ifndef STABRK3_H_
