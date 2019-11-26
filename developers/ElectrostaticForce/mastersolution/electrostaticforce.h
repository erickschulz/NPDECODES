/**
 * @file newproblem.h
 * @brief NPDE homework NewProblem code
 * @author Oliver Rietmann, Erick Schulz
 * @date 01.01.2020
 * @copyright Developed at ETH Zurich
 */

#include <iostream>
// Eigen includes
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
// Lehrfem++ includes
#include <lf/assemble/assemble.h>
#include <lf/geometry/geometry.h>
#include <lf/io/io.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/utils/utils.h>
#include <lf/uscalfe/uscalfe.h>

namespace ElectrostaticForce {

/** @brief Compute the electrostatic force
 *       F(u) = (1/2) * int_{\Gamma_1} <grad(u),n> * grad(u) dS
 * from the potential directly using the trapezoidal rule. Recall that
 * for a partition x_0 < x_1 < ... < x_N of a curve C. The trapezoidal rule
 * reads which approximates
 *                         int_C f(x) dS
 * reads
 *      ( dx/2 ) * ( f(x_0) + 2*f(x_1) + ... + 2*f(x_N-1) + f(X_N) ).
 * Hence, when the curve is closed, that is when the endpoints x_0 = x_N are
 * the same, the quadrature rule simply reads
 *                      dx * sum_{i=0}^N f(x_i).
 * @return two dimensional force vector */
Eigen::Vector2d computeExactForce();

/** @Brief This function enforces Dirichlet zero boundary conditions on the
 * Galerkin matrices. It transforms every columns and rows associated to a
 * global index belonging to a degree of freedom lying on the boundary to zero
 * entries but the diagonal one which is set to 1.0
 * @param selectvals is the predicate identifying the boundary indices of the
 * rows and columns that are to be dropped */
template <typename SCALAR, typename SELECTOR>
void dropMatrixRowsAndColumns(SELECTOR &&selectvals,
                              lf::assemble::COOMatrix<SCALAR> &A) {
  const lf::assemble::size_type N(A.cols());
  LF_ASSERT_MSG(A.rows() == N, "Matrix must be square!");
  A.setZero(
      [&selectvals](lf::assemble::gdof_idx_t i, lf::assemble::gdof_idx_t j) {
        return (selectvals(i) || selectvals(j));
      });
  for (lf::assemble::gdof_idx_t dofnum = 0; dofnum < N; ++dofnum) {
    const auto selval{selectvals(dofnum)};
    if (selval) {
      A.AddToEntry(dofnum, dofnum, 1.0);
    }
  }
}

Eigen::VectorXd solvePoissonBVP(
    const std::shared_ptr<lf::uscalfe::FeSpaceLagrangeO1<double>> &fe_space_p);

}  // namespace ElectrostaticForce
