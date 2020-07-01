/**
 * @file engquistoshernumericalflux.h
 * @brief NPDE homework "EngquistOsherNumericalFlux" code
 * @author Oliver Rietmann
 * @date 23.04.2019
 * @copyright Developed at ETH Zurich
 */

#include <Eigen/Core>

namespace EngquistOsherNumericalFlux {

/**
 * @brief Engquist-Osher numerical flux F_EO(v, w)
 * for the flux function f(u) = cosh(u).
 *
 * @param v left state: real number satisfying v <= w
 * @param w right state: real number satisfying v <= w
 * @return F_EO(v, w)
 */
double EngquistOsherNumFlux(double v, double w);

/**
 * @brief Computes the solution u(x, T) for inital
 * data u0 on [a, b] (assuming continuous, constant
 * continuation to the reals) under the
 * Engquist-Osher numerical flux from above.
 *
 * @param a left bound of spacial interval, a < b
 * @param b right bound of spacial interval, a < b
 * @return Vector of same size as u0, containing
 * the (approximate) values of u(x, T) on the same
 * spacial grid as u0.
 */
Eigen::VectorXd solveCP(double a, double b, Eigen::VectorXd u0, double T);

} // namespace EngquistOsherNumericalFlux
