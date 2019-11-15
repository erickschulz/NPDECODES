/**
 * @file burgersequation.h
 * @brief NPDE homework BurgersEquation code
 * @author Oliver Rietmann
 * @date 15.04.2019
 * @copyright Developed at ETH Zurich
 */

#include <Eigen/Core>

namespace BurgersEquation {

const static Eigen::IOFormat CSVFormat(Eigen::FullPrecision,
                                       Eigen::DontAlignCols, ", ", "\n");

/**
 * @brief Computes the solution u(x, t) on [-1, 4]x[0, T] of
 * Burgers equation using the Godunov numerical flux function
 * for discretization in space and the explicit Euler scheme
 * for discretization in time, given the inital data
 * u(x, 0) = sin(PI * x) for x in [0, 1] and zero otherwise.
 *
 * @param T final time up to which we approximate the solution.
 * @param N number of sub intervals of [-1, 4] of the discretization.
 * @return nodal values of the approximate solution at time T,
 * where the timesteps are of same size as the spacial mesh-width.
 */
Eigen::VectorXd solveBurgersGodunov(double T, unsigned int N);

/**
 * @brief Computes the error of solveBurgersGodunov(T, N) with the
 * spacial mesh-width (= timestep size) for T = 0.3 and T = 3.0.
 * The reference solution is computed with solveBurgersGodunov(T, N)
 * for N = 3200, in the the time-dependent, descrete L1-norm.
 *
 * @return The first row contains the values
 * h = {0.1, 0.05, 0.025, 0.00125}. The second and third row contain
 * the error for T = 0.3 and T = 3.0 respectively.
 */
Eigen::Matrix<double, 3, 4> numexpBurgersGodunov();

}  // namespace BurgersEquation
