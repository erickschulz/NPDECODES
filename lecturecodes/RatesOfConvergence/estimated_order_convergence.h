#ifndef RATE_CVG_H
#define RATE_CVG_H

/** @file
 * @brief functions for estimating order of convergence
 * @author Ralf Hiptmair
 * @date March 2019
 */

#include <string>

#include <Eigen/Core>

namespace rate_of_convergence {

Eigen::Vector2d linearFit(const Eigen::VectorXd& x, const Eigen::VectorXd& y);

double eoc(const Eigen::VectorXd& N, const Eigen::VectorXd& err,
           unsigned fromindex = 0, std::string filename = "conv.eps");

void plotConvergence(const Eigen::VectorXd& N, const Eigen::VectorXd& err,
                     const Eigen::VectorXd& fit, const std::string& label,
                     const std::string& filename);

}  // namespace rate_of_convergence

#endif
