///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// (C) 2016 SAM, D-MATH
/// Author(s): Julien Gacon
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

# ifndef POLYFIT_HPP
# define POLYFIT_HPP

# include <Eigen/Dense>
# include <Eigen/QR>

/* SAM_LISTING_BEGIN_0 */
// Solver for polynomial linear least squares data fitting problem
// data points passed in t and y, 'order' = degree of polynomial 
Eigen::VectorXd polyfit(const Eigen::VectorXd& t, const Eigen::VectorXd& y, const unsigned& order) {
  // A = [1 t_1 t_1^2 ... ]
  //     [ ...        ... ]
  //     [1 t_n t_n^2 ... ]
  Eigen::MatrixXd A = Eigen::MatrixXd::Ones(t.size(),order+1);
  for (unsigned j = 1; j <= order; ++j) {
    A.col(j) = A.col(j - 1).cwiseProduct(t);
  }
  Eigen::VectorXd coeffs = A.householderQr().solve(y);

  return coeffs.reverse();
}
/* SAM_LISTING_END_0 */

# endif
