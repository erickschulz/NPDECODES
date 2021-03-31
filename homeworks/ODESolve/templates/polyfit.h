///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Julien Gacon
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

#ifndef POLYFIT_HPP
#define POLYFIT_HPP

#include <Eigen/Core>

namespace ODESolve {

/* SAM_LISTING_BEGIN_0 */
// Solver for polynomial linear least squares data fitting problem
// data points passed in t and y, 'order' = degree of polynomial
Eigen::VectorXd polyfit(const Eigen::VectorXd& t, const Eigen::VectorXd& y,
                        const unsigned& order);
/* SAM_LISTING_END_0 */

}  // namespace ODESolve

#endif
