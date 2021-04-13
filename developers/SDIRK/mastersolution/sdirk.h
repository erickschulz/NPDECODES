#ifndef SDIRK_H_
#define SDIRK_H_

/**
 * @file sdirk.h
 * @brief NPDE homework SDIRK code
 * @author Unknown, Oliver Rietmann
 * @date 31.03.2021
 * @copyright Developed at ETH Zurich
 */

#include <Eigen/Core>
#include <vector>

namespace SDIRK {

/**
 * @brief  Compute one step of  the IVP y'' + y' + y = 0 using a SDIRK method
 * @param z0 initial state z0 = [y(0), y'(0)]
 * @param h size of the step
 * @param gamma parameter of the scheme
 * @return next step z1
 */
Eigen::Vector2d SdirkStep(const Eigen::Vector2d &z0, double h, double gamma);

/**
 * @brief Solve autonomous IVP y'' + y' + y = 0, [y(0), y'(0)] = z0 using SDIRK
 * method with N equidistant steps
 * @param z0  initial data z0 = [y(0), y'(0)]
 * @param N number of equidistant steps
 * @param T final time
 * @param gamma parameter of the scheme
 * @return vector containing each step z_k (y and y')
 */
std::vector<Eigen::Vector2d> SdirkSolve(const Eigen::Vector2d &z0,
                                        unsigned int N, double T, double gamma);

double CvgSDIRK();

}  // namespace SDIRK

#endif  // #define SDIRK_H_
