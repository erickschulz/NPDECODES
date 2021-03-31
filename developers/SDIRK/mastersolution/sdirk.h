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

//! \brief One step of autonomous IVP y'' + y' + y = 0, [y(0), y'(0)] = z0 using
//! SDIRK method Use SDIRK method for first order ode z' = f(z). Steps of size
//! h.
//! \tparam Eigen::VectorXd type of solution space y and initial data y0
//! \param[in] z0 initial data z(0)
//! \param[in] h size of the step
//! \param[in] gamma parameter
//! \return next step z1
Eigen::Vector2d sdirkStep(const Eigen::Vector2d &z0, double h, double gamma);

//! \brief Solve autonomous IVP y'' + y' + y = 0, [y(0), y'(0)] = z0 using SDIRK
//! method. Use SDIRK method for first order ode z' = f(z), with N equidistant
//! steps.
//! \tparam Eigen::VectorXd type of solution space z = [y,y']! and initial
//! data z0 = [y(0), y'(0)] \param[in] z0 initial data z(0)
//! \param[in] N number of equidistant steps
//! \param[in] T final time of simulation
//! \param[in] gamma parameter
//! \return vector containing each step of z_k (y and y')
std::vector<Eigen::Vector2d> sdirkSolve(const Eigen::Vector2d &z0,
                                        unsigned int N, double T, double gamma);

double cvgSDIRK();

}  // namespace SDIRK

#endif  // #define SDIRK_H_
